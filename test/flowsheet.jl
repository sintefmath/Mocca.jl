using Test
using Mocca
using Jutul

@testset "Device laws and profiles" begin
    R = Mocca.GAS_CONSTANT
    st = (InletPressure = [2e5], OutletPressure = [1e5], InletTemperature = [300.0], OutletTemperature = [350.0])

    # Linear valve uses the upstream concentration
    v = Mocca.LinearValve(1e-6; opening = 0.5)
    @test Mocca.device_flow(v, st, 0.0) ≈ 0.5 * 1e-6 * 1e5 * 2e5 / (R * 300.0)
    st_rev = (InletPressure = [1e5], OutletPressure = [2e5], InletTemperature = [300.0], OutletTemperature = [350.0])
    @test Mocca.device_flow(v, st_rev, 0.0) ≈ -0.5 * 1e-6 * 1e5 * 2e5 / (R * 350.0)

    # Flow controller: volumetric rate at outlet pressure and inlet temperature
    @test Mocca.device_flow(Mocca.VolumetricFlow(0.01), st, 0.0) ≈ 0.01 * 1e5 / (R * 300.0)
    @test Mocca.device_flow(Mocca.Closed(), st, 0.0) == 0

    ramp = Mocca.ExponentialRamp(1e4, 1e5, 0.5)
    @test Mocca.evaluate_profile(ramp, 0.0) ≈ 1e4
    @test Mocca.evaluate_profile(ramp, 100.0) ≈ 1e5
    shifted = Mocca.shift_profile(ramp, 10.0)
    @test Mocca.evaluate_profile(shifted, 10.0) ≈ 1e4
    @test Mocca.evaluate_profile(3.0, 42.0) == 3.0

    @test Mocca.light_product_composition([0.15, 0.85]) ≈ [1e-10, 1 - 1e-10]
end

@testset "Flowsheet validation" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    bed = Mocca.setup_process_model(Mocca.FixedBed(constants), constants; ncells = 5)
    names = Mocca.component_names(bed.system)
    fs = Mocca.Flowsheet()
    Mocca.add_unit!(fs, :Bed, bed)
    Mocca.add_unit!(fs, :V1, Mocca.setup_device_model(names))
    Mocca.add_unit!(fs, :V2, Mocca.setup_device_model(names))
    @test_throws ErrorException Mocca.add_unit!(fs, :Bed, bed)
    @test_throws ErrorException Mocca.connect!(fs, (:V1, :outlet) => (:V2, :inlet))
    @test_throws ErrorException Mocca.connect!(fs, (:Bed, :middle) => (:V1, :outlet))

    # Either order is accepted
    Mocca.connect!(fs, (:V1, :outlet) => (:Bed, :bottom))
    @test_throws ErrorException Mocca.connect!(fs, (:Bed, :top) => (:V1, :outlet))
    Mocca.connect!(fs, (:Bed, :top) => (:V2, :inlet))
    @test Set(Mocca.boundary_ports(fs)) == Set([(:V1, :inlet), (:V2, :outlet)])

    # Unconnected ports need a boundary condition
    @test_throws ErrorException Mocca.setup_flowsheet_model(fs)
    bc = Mocca.PortCondition(pressure = 1e5, temperature = 298.15, composition = [0.15, 0.85])
    Mocca.set_boundary!(fs, :V1, :inlet, bc)
    Mocca.set_boundary!(fs, :V2, :outlet, bc)
    mm = Mocca.setup_flowsheet_model(fs)
    @test mm isa Jutul.MultiModel
    # Two cross terms per connection plus one for each energy balance
    @test length(mm.cross_terms) == 6

    # A stage cannot set a condition on a connected port
    bad = [Mocca.Stage("bad", 1.0; V1 = (outlet = bc,))]
    @test_throws ErrorException Mocca.setup_schedule(fs, mm, bad)
end

function _vsa_legacy(constants, init; ncells, num_cycles)
    bed = Mocca.setup_process_model(Mocca.FixedBed(constants), constants; ncells = ncells)
    state0 = Mocca.setup_process_state(bed; init...)
    prm = Mocca.setup_process_parameters(bed)
    bcs = Mocca.setup_boundary_conditions(constants, ["pressurisation", "adsorption", "blowdown", "evacuation"])
    forces, dt = Mocca.setup_forces(bed, [15.0, 15.0, 30.0, 40.0], bcs; num_cycles = num_cycles, max_dt = 1.0)
    states, timesteps = Mocca.simulate_process(Mocca.MoccaCase(bed, dt, forces; state0 = state0, parameters = prm);
        info_level = -1, output_substates = true)
    return (states, timesteps)
end

function _vsa_flowsheet(constants, init; ncells, num_cycles)
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = ncells)
    mm = Mocca.setup_flowsheet_model(fs)
    state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.setup_process_state(fs[:Bed]; init...))
    prm = Mocca.setup_flowsheet_parameters(fs)
    forces, dt = Mocca.setup_schedule(fs, mm, stages; num_cycles = num_cycles, max_dt = 1.0)
    states, timesteps = Mocca.simulate_process(Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = prm);
        info_level = -1, output_substates = true)
    return (fs, states, timesteps, state0, prm)
end

# State at the end of the step finishing at time t
_state_at(states, timesteps, t) = states[findfirst(≈(t), cumsum(timesteps))]

@testset "Legacy inlet boundary conditions carry the upstream composition" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    model = Mocca.setup_process_model(Mocca.FixedBed(constants), constants; ncells = 10)
    pressurisation, adsorption = Mocca.setup_boundary_conditions(constants, ["pressurisation", "adsorption"])
    y_cell = [0.02, 0.98]
    cell_state(P) = (Pressure = fill(P, 10), Temperature = fill(298.15, 10), y = repeat(y_cell, 1, 10))
    # Negative flux is flow into the column. At the start of pressurisation the
    # boundary is at the low pressure.
    for (bc, P, into_column, y_expected) in (
            (adsorption, constants.p_high, true, constants.y_feed),
            (pressurisation, 0.5 * constants.p_low, true, constants.y_feed),
            (pressurisation, 2.0 * constants.p_low, false, y_cell),
        )
        flux = Mocca.mass_flux_left(cell_state(P), model, 0.0, bc)
        @test (sum(flux) < 0) == into_column
        @test flux ./ sum(flux) ≈ y_expected
    end
end

@testset "Flowsheet VSA matches the legacy boundary conditions" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    # Start from an evacuated bed, as in a running cycle
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    ncells = 30
    legacy, legacy_dt = _vsa_legacy(constants, init; ncells = ncells, num_cycles = 1)
    fs, states, timesteps, state0, prm = _vsa_flowsheet(constants, init; ncells = ncells, num_cycles = 1)
    @test sum(timesteps) ≈ sum(legacy_dt) ≈ 100.0

    for t in (15.0, 30.0, 60.0, 100.0)  # end of each stage
        a = _state_at(legacy, legacy_dt, t)
        b = _state_at(states, timesteps, t)[:Bed]
        @test isapprox(b[:Pressure], a[:Pressure], rtol = 1e-3)
        @test isapprox(maximum(b[:Temperature]), maximum(a[:Temperature]), atol = 0.1)
        @test isapprox(sum(b[:AdsorbedConcentration][1, :]), sum(a[:AdsorbedConcentration][1, :]), rtol = 2e-3)
        @test isapprox(sum(b[:AdsorbedConcentration][2, :]), sum(a[:AdsorbedConcentration][2, :]), rtol = 2e-3)
    end

    # Total moles are conserved: what the devices move in and out equals the
    # change of gas and adsorbed inventory in the bed
    R = Mocca.GAS_CONSTANT
    bp = prm[:Bed]
    inventory(s) = sum(s[:Pressure] ./ (R .* s[:Temperature]) .* bp[:FluidVolume]) +
        sum(sum(s[:AdsorbedConcentration], dims = 1)' .* bp[:SolidVolume])
    net_in = 0.0
    for (s, dt) in zip(states, timesteps)
        net_in += (s[:V_feed][:MolarFlow][1] - s[:V_product][:MolarFlow][1] - s[:V_vacuum][:MolarFlow][1]) * dt
    end
    ΔM = inventory(states[end][:Bed]) - inventory(state0[:Bed])
    throughput = sum(abs(s[:V_feed][:MolarFlow][1]) * dt for (s, dt) in zip(states, timesteps))
    @test abs(net_in - ΔM) < 1e-8 * throughput
end

@testset "Flowsheet plotting" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    legacy, legacy_dt = _vsa_legacy(constants, init; ncells = 10, num_cycles = 1)
    fs, states, timesteps, = _vsa_flowsheet(constants, init; ncells = 10, num_cycles = 1)
    model = Mocca.setup_flowsheet_model(fs)
    _, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 10)
    bed = fs[:Bed]
    bed_states = [s[:Bed] for s in states]
    runs = ["Flowsheet" => (bed_states, timesteps), "Boundary conditions" => (legacy, legacy_dt)]
    @test Mocca.plot_cell_comparison(runs, bed, 10) isa Mocca.Figure
    @test Mocca.plot_state_comparison(["Flowsheet" => bed_states[end], "Boundary conditions" => legacy[end]], bed) isa Mocca.Figure
    @test Mocca.plot_flowsheet_streams(states, model, timesteps; stages = stages) isa Mocca.Figure
    @test Mocca.plot_flowsheet_streams(states, model, timesteps) isa Mocca.Figure
end

@testset "Stage time steps" begin
    # Without first_dt, equal steps no longer than max_dt, as before
    @test Mocca._stage_steps(15.0, 2.0, 2.0) ≈ fill(15.0 / 8, 8)
    s = Mocca._stage_steps(20.0, 1.0, 0.01)
    @test sum(s) ≈ 20.0
    @test s[1:3] ≈ [0.01, 0.02, 0.04]
    @test maximum(s) <= 1.0 + 1e-12
    @test Mocca._stage_steps(0.005, 1.0, 0.01) ≈ [0.005]
end

@testset "Parallel stage sequences" begin
    cond(P) = Mocca.PortCondition(pressure = P, temperature = 298.15, composition = [0.5, 0.5])
    ramp = Mocca.ExponentialRamp(1e5, 2e4, 0.5)
    a = [Mocca.Stage("feed", 10.0; V_a = (law = Mocca.Closed(),)),
         Mocca.Stage("vent", 20.0; V_a = (law = Mocca.LinearValve(1e-6), outlet = cond(ramp)))]
    b = [Mocca.Stage("feed", 10.0; V_b = (law = Mocca.Closed(),)),
         Mocca.Stage("vent", 5.0; V_b = (law = Mocca.Closed(),)),
         Mocca.Stage("evacuate", 25.0; V_b = (law = Mocca.Closed(),))]
    stages = Mocca.parallel_stages("A" => a, "B" => b)
    @test [st.duration for st in stages] ≈ [10.0, 5.0, 15.0, 10.0]
    @test [st.name for st in stages] == ["A feed, B feed", "A vent, B vent", "A vent, B evacuate", "A idle, B evacuate"]
    @test Set(keys(stages[3].settings)) == Set([:V_a, :V_b])
    @test isempty(setdiff(keys(stages[4].settings), [:V_b]))
    # A's vent stage is split after 5 s; its ramp keeps its own start time
    split_ramp = Mocca.shift_profile(stages[3].settings[:V_a].outlet, 15.0).pressure
    whole_ramp = Mocca.shift_profile(a[2].settings[:V_a].outlet, 10.0).pressure
    for t in (15.0, 22.0, 30.0)
        @test Mocca.evaluate_profile(split_ramp, t) ≈ Mocca.evaluate_profile(whole_ramp, t)
    end
    @test_throws ErrorException Mocca.parallel_stages("A" => a, "A again" => a)
    @test_throws ErrorException Mocca.parallel_stages("A" => a, "B" => b; cycle_time = 20.0)
end

@testset "Wet flue gas flowsheets" begin
    function run_one_cycle(fs, stages, state0)
        model = Mocca.setup_flowsheet_model(fs)
        prm = Mocca.setup_flowsheet_parameters(fs)
        forces, dt = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0)
        s0 = Mocca.setup_flowsheet_state(fs; state0...)
        case = Mocca.MoccaCase(model, dt, forces; state0 = s0, parameters = prm)
        states, ts = Mocca.simulate_process(case; info_level = -1, output_substates = true)
        @test sum(ts) ≈ sum(dt)
        return (states, ts, s0, prm)
    end
    # Each component is conserved: what the boundary devices move in and out
    # equals the change of inventory in the beds
    function check_conservation(fs, states, ts, s0, prm, beds, inflow, outflow)
        inventory(s) = sum(Mocca.component_inventory(fs[b], s[b], prm[b]) for b in beds)
        net_in = sum(Mocca.stream_totals(states, ts, d) for d in inflow) .- sum(Mocca.stream_totals(states, ts, d) for d in outflow)
        Δ = inventory(states[end]) .- inventory(s0)
        fed = Mocca.stream_totals(states, ts, :V_feed)
        @test all(abs.(net_in .- Δ) .< 1e-6 * sum(fed))
    end

    fs, stages = Mocca.lpp_vsa_flowsheet(ncells = 10)
    @test [st.name for st in stages] == ["LPP", "adsorption", "blowdown", "evacuation"]
    states, ts, s0, prm = run_one_cycle(fs, stages, (; Bed = Mocca.wet_flue_gas_initial_state(fs[:Bed], 0.03e5)))
    check_conservation(fs, states, ts, s0, prm, (:Bed,), (:V_feed, :V_lpp), (:V_product, :V_vacuum))

    # Cycling until the CO2 balance closes: with a loose tolerance it stops as
    # soon as min_cycles is reached, with a tight one it runs to max_cycles
    model = Mocca.setup_flowsheet_model(fs)
    forces, dt = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 2.0, first_dt = 0.1)
    s0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.wet_flue_gas_initial_state(fs[:Bed], 0.03e5))
    kw = (inflow = (:V_feed, :V_lpp), outflow = (:V_product, :V_vacuum), consecutive = 1, min_cycles = 2)
    loose = Mocca.simulate_until_mass_balance(model, s0, prm, forces, dt; kw..., tol = 1.0, max_cycles = 3)
    @test loose.converged && loose.cycles == 2
    tight = Mocca.simulate_until_mass_balance(model, s0, prm, forces, dt; kw..., tol = 1e-12, max_cycles = 2)
    @test !tight.converged && tight.cycles == 2 && length(tight.history) == 2
    @test tight.history[2] < tight.history[1]

    fs, stages = Mocca.dual_adsorbent_vsa_flowsheet(ncells = 10)
    @test sum(st -> st.duration, stages) ≈ 20.0 + 46.20 + 56.30 + 101.20
    states, ts, s0, prm = run_one_cycle(fs, stages, (;
        SilicaGel = Mocca.wet_flue_gas_initial_state(fs[:SilicaGel], 0.30e5),
        Zeolite = Mocca.wet_flue_gas_initial_state(fs[:Zeolite], 101325.0)))
    check_conservation(fs, states, ts, s0, prm, (:SilicaGel, :Zeolite), (:V_feed, :V_lpp), (:V_waste, :V_product, :V_vacuum))
    # The silica gel bed holds back the water
    n_feed = Mocca.stream_totals(states, ts, :V_feed)
    n_link = Mocca.stream_totals(states, ts, :V_link)
    @test n_link[3] < 1e-3 * n_feed[3]
end

@testset "Two-bed VSA with pressure equalisation" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    fs, stages = Mocca.two_bed_vsa_flowsheet(constants; ncells = 20)
    @test length(stages) == 6
    mm = Mocca.setup_flowsheet_model(fs)
    init(P) = Mocca.setup_process_state(fs[:A]; Pressure = P, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    state0 = Mocca.setup_flowsheet_state(fs; A = init(constants.p_high / 2), B = init(constants.p_low))
    prm = Mocca.setup_flowsheet_parameters(fs)
    forces, dt = Mocca.setup_schedule(fs, mm, stages; num_cycles = 1, max_dt = 1.0)
    states, timesteps = Mocca.simulate_process(Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = prm);
        info_level = -1, output_substates = true)
    @test sum(timesteps) ≈ 80.0

    P_mean(s, b) = sum(s[b][:Pressure]) / length(s[b][:Pressure])
    before = _state_at(states, timesteps, 30.0)
    after = _state_at(states, timesteps, 40.0)
    # A gives gas to B through the equalisation valve, closing most of the gap
    @test after[:V_eq][:MolarFlow][1] > 0
    gap_before = P_mean(before, :A) - P_mean(before, :B)
    gap_after = P_mean(after, :A) - P_mean(after, :B)
    @test gap_before > 0.5 * constants.p_high
    @test abs(gap_after) < 0.2 * gap_before
    # In the second equalisation B gives gas back to A
    @test _state_at(states, timesteps, 80.0)[:V_eq][:MolarFlow][1] < 0

    # Moles are conserved across both beds; the equalisation flow is internal
    R = Mocca.GAS_CONSTANT
    inventory(s, b) = sum(s[b][:Pressure] ./ (R .* s[b][:Temperature]) .* prm[b][:FluidVolume]) +
        sum(sum(s[b][:AdsorbedConcentration], dims = 1)' .* prm[b][:SolidVolume])
    net_in = 0.0
    throughput = 0.0
    for (s, h) in zip(states, timesteps)
        for b in (:A, :B)
            F_feed = s[Symbol(:V_feed_, b)][:MolarFlow][1]
            net_in += (F_feed - s[Symbol(:V_product_, b)][:MolarFlow][1] - s[Symbol(:V_vacuum_, b)][:MolarFlow][1]) * h
            throughput += abs(F_feed) * h
        end
    end
    ΔM = sum(inventory(states[end], b) - inventory(state0, b) for b in (:A, :B))
    @test abs(net_in - ΔM) < 1e-8 * throughput
end
