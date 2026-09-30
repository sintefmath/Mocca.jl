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

function _vsa_flowsheet(constants, init; ncells, num_cycles, legacy_inlet = false)
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = ncells)
    mm = Mocca.setup_flowsheet_model(fs; legacy_inlet = legacy_inlet)
    state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.setup_process_state(fs[:Bed]; init...))
    prm = Mocca.setup_flowsheet_parameters(fs)
    forces, dt = Mocca.setup_schedule(fs, mm, stages; num_cycles = num_cycles, max_dt = 1.0)
    states, timesteps = Mocca.simulate_process(Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = prm);
        info_level = -1, output_substates = true)
    return (fs, states, timesteps, state0, prm)
end

# State at the end of the step finishing at time t
_state_at(states, timesteps, t) = states[findfirst(≈(t), cumsum(timesteps))]

@testset "Flowsheet VSA matches the legacy boundary conditions" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    # Start from an evacuated bed, as in a running cycle
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    ncells = 30
    legacy, legacy_dt = _vsa_legacy(constants, init; ncells = ncells, num_cycles = 1)
    fs, states, timesteps, state0, prm = _vsa_flowsheet(constants, init; ncells = ncells, num_cycles = 1, legacy_inlet = true)
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
