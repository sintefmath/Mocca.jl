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
    bed = haghpanah_bed(ncells = 5)
    names = Mocca.component_names(bed.system)
    fs = Mocca.UnitNetwork()
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
    @test Set(Mocca.boundary_ports(fs)) == Set([Mocca.PortRef(:V1, :inlet), Mocca.PortRef(:V2, :outlet)])

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

# A connection type defined outside Mocca, joining two unit ports directly and
# counting the cross terms it is asked for
struct _TestTie <: Mocca.AbstractConnection
    a::Mocca.PortRef
    b::Mocca.PortRef
    calls::Base.RefValue{Int}
end
Mocca.connection_ports(c::_TestTie) = (c.a, c.b)
function Mocca.check_connection(net, c::_TestTie)
    all(r -> haskey(net.units, r.name), (c.a, c.b)) || error("unknown unit")
    any(x -> x isa _TestTie, net.connections) && error("one tie only")
end
Mocca.add_cross_terms!(mm, net, c::_TestTie) = (c.calls[] += 1; mm)

@testset "Connections are an open type" begin
    names = ["CO2", "N2"]
    net = Mocca.UnitNetwork()
    for b in (:A, :B)
        Mocca.add_unit!(net, b, haghpanah_bed(ncells = 4))
    end
    Mocca.add_unit!(net, :V, Mocca.setup_device_model(names))
    Mocca.connect!(net, (:A, :top) => (:V, :inlet))
    Mocca.connect!(net, Mocca.StreamConnection(Mocca.PortRef(:B, :bottom), Mocca.PortRef(:V, :outlet)))
    calls = Ref(0)
    @test_throws ErrorException Mocca.connect!(net, _TestTie(Mocca.PortRef(:A, :bottom), Mocca.PortRef(:C, :top), calls))
    Mocca.connect!(net, _TestTie(Mocca.PortRef(:A, :bottom), Mocca.PortRef(:B, :top), calls))
    @test_throws ErrorException Mocca.connect!(net, _TestTie(Mocca.PortRef(:A, :bottom), Mocca.PortRef(:B, :top), calls))
    @test length(Mocca.connections(net)) == 3 && length(Mocca.stream_connections(net)) == 2
    @test Mocca.is_connected(net, :A, :bottom) && Mocca.is_connected(net, :V, :outlet)
    @test Mocca.attachments(net, :V) == Mocca.stream_connections(net)
    mm = Mocca.setup_flowsheet_model(net)
    # Each connection adds its own cross terms: the tie none, each stream
    # three (port state, mass and energy)
    @test calls[] == 1
    @test length(mm.cross_terms) == 2 * 3
    @test occursin("UnitNetwork with 2 units and 1 device", sprint(show, MIME"text/plain"(), net))
end

@testset "Streams carry the enthalpy of the upstream port" begin
    bed = haghpanah_bed(ncells = 3)
    sys = bed.system
    eq_m, eq_e = bed.equations[:mass_conservation], bed.equations[:energy_column]
    device(F, y_in, T_in, y_out, T_out; P_in = 1e5) = (MolarFlow = [F],
        InletComposition = reshape(y_in, :, 1), InletTemperature = [T_in], InletPressure = [P_in],
        OutletComposition = reshape(y_out, :, 1), OutletTemperature = [T_out], OutletPressure = [0.5e5])
    cp(y, T = 300.0) = Mocca.molar_heat_capacity(sys.thermodynamics, y, T)
    cell(T, y) = (Temperature = fill(T, 3), y = repeat(y, 1, 3), C_pg = fill(cp(y, T), 3))
    function term(ct, eq, target, dev; model = bed)
        out = zeros(eq isa typeof(eq_m) ? 2 : 1)
        Jutul.update_cross_term_in_entity!(out, 1, target, target, dev, dev, model, nothing, ct, eq, 1.0)
        return out
    end
    # A link from the top of bed A (cell 3) to the bottom of bed B (cell 1),
    # whose cells differ in temperature and composition
    yA, TA, yB, TB = [0.3, 0.7], 320.0, [0.05, 0.95], 300.0
    A, B = cell(TA, yA), cell(TB, yB)
    for F in (0.02, -0.03)
        dev = device(F, yA, TA, yB, TB)
        y_up, T_up = F > 0 ? (yA, TA) : (yB, TB)
        # What leaves one bed enters the other, component by component
        @test term(Mocca.StreamCT(3, 1), eq_m, A, dev) ≈ -term(Mocca.StreamCT(1, 2), eq_m, B, dev)
        # The bed it leaves takes C_pg·T, which leaves its temperature as it
        # is, and the bed it enters gets C_pg·T of its own gas and the heat to
        # bring the stream to its temperature
        (Tl, cl, Te, ce, out, in) = F > 0 ? (TA, A, TB, B, Mocca.StreamCT(3, 1), Mocca.StreamCT(1, 2)) : (TB, B, TA, A, Mocca.StreamCT(1, 2), Mocca.StreamCT(3, 1))
        @test term(out, eq_e, cl, dev)[1] ≈ abs(F) * cl.C_pg[1] * Tl
        @test term(in, eq_e, ce, dev)[1] ≈ -abs(F) * (ce.C_pg[1] * Te + cp(y_up) * (T_up - Te))
    end
    # A stream from a boundary brings the same, whatever its pressure; the
    # molar heat capacities of the constants are nearly equal (30.7 J/(mol·K)),
    # so nearly its own c_p·T
    feed = [0.15, 0.85]
    for (T, y, P) in ((300.0, [1e-10, 1 - 1e-10], 1e5), (350.0, [0.5, 0.5], 0.1e5))
        c = cell(T, y)
        e = term(Mocca.StreamCT(1, 2), eq_e, c, device(0.1, feed, 298.15, feed, 298.15; P_in = P))
        @test e ≈ [-0.1 * (c.C_pg[1] * T + cp(feed) * (298.15 - T))]
        @test e ≈ [-0.1 * cp(feed) * 298.15] rtol = 1e-4
    end

    # With heat capacities that change with temperature the stream brings the
    # integral of the heat capacity, not the difference of c_p·T
    k = Mocca.convert_constants(haghpanah_constants())
    gas = Gas(["CO2", "N2"]; molar_masses = k.gas.molar_masses, thermodynamics = :shomate)
    hot = Mocca.fixed_bed_model(gas, k.feed, k.adsorbent, k.column; ncells = 3, dispersion_velocity = k.dispersion_velocity)
    th = hot.system.thermodynamics
    Tc, Tin = 360.0, 298.15
    st = (Temperature = fill(Tc, 3), y = repeat(feed, 1, 3), C_pg = fill(Mocca.molar_heat_capacity(th, feed, Tc), 3))
    e = term(Mocca.StreamCT(1, 2), hot.equations[:energy_column], st, device(0.1, feed, Tin, feed, Tin); model = hot)
    # Trapezoidal integral of the heat capacity from the stream's temperature to the cell's
    Ts = range(Tin, Tc, length = 2001)
    ∫cp = sum(Mocca.molar_heat_capacity(th, feed, T) for T in Ts[2:end-1]) * step(Ts) +
        (Mocca.molar_heat_capacity(th, feed, Tin) + Mocca.molar_heat_capacity(th, feed, Tc)) * step(Ts) / 2
    @test e[1] ≈ -0.1 * (st.C_pg[1] * Tc - ∫cp) rtol = 1e-8
end

# State at the end of the step finishing at time t
_state_at(states, timesteps, t) = states[findfirst(≈(t), cumsum(timesteps))]

@testset "Flowsheet VSA conserves moles" begin
    # From an evacuated bed, as in a running cycle
    p = haghpanah_vsa(; ncells = 30)
    state0, prm = Mocca.initial_state(p), p.parameters
    r = simulate(p; state0 = state0)
    states, timesteps = r.states, r.timesteps
    @test sum(timesteps) ≈ 100.0
    # The stages end where they should, with the pressures they ramp to
    @test maximum(_state_at(states, timesteps, 30.0)[:Bed][:Pressure]) ≈ 1e5 rtol = 0.05
    @test minimum(_state_at(states, timesteps, 100.0)[:Bed][:Pressure]) < 0.2e5

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
    p = haghpanah_vsa(; ncells = 10)
    q = haghpanah_vsa(haghpanah_constants(v_feed = 0.3); ncells = 10)
    r, rq = simulate(p), simulate(q)
    bed = p.flowsheet[:Bed]
    runs = ["v_feed 0.37 m/s" => (unit_states(r.states), r.timesteps), "v_feed 0.3 m/s" => (unit_states(rq.states), rq.timesteps)]
    @test Mocca.plot_cell_comparison(runs, bed, 10) isa Mocca.Plotting.Figure
    @test Mocca.plot_state_comparison(["v_feed 0.37 m/s" => r.states[end][:Bed], "v_feed 0.3 m/s" => rq.states[end][:Bed]], bed) isa Mocca.Plotting.Figure
    @test Mocca.plot_flowsheet_streams(r.states, p.model, r.timesteps; stages = Mocca.stages(p)) isa Mocca.Plotting.Figure
    @test Mocca.plot_flowsheet_streams(r.states, p.model, r.timesteps) isa Mocca.Plotting.Figure
end

@testset "Stage time steps" begin
    # Without first_dt, equal steps no longer than max_dt, as before
    @test Mocca.Operation._stage_steps(15.0, 2.0, 2.0) ≈ fill(15.0 / 8, 8)
    s = Mocca.Operation._stage_steps(20.0, 1.0, 0.01)
    @test sum(s) ≈ 20.0
    @test s[1:3] ≈ [0.01, 0.02, 0.04]
    @test maximum(s) <= 1.0 + 1e-12
    @test Mocca.Operation._stage_steps(0.005, 1.0, 0.01) ≈ [0.005]
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

@testset "Wet flue gas processes" begin
    ts = TimeStepping(max_dt = 1.0, first_dt = 0.01)
    function run_one_cycle(p)
        s0 = Mocca.initial_state(p)
        r = simulate(p; state0 = s0)
        return (r.states, r.timesteps, s0, p.parameters)
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

    p = lpp_vsa_process(ncells = 10, timestepping = ts)
    @test [st.name for st in Mocca.stages(p)] == ["LPP", "adsorption", "blowdown", "evacuation"]
    states, steps, s0, prm = run_one_cycle(p)
    check_conservation(p.flowsheet, states, steps, s0, prm, (:Bed,), (:V_feed, :V_lpp), (:V_product, :V_vacuum))

    # Cycling until the CO2 balance closes: with a loose tolerance it stops as
    # soon as min_cycles is reached, with a tight one it runs to max_cycles
    forces, dt = Mocca.setup_schedule(p.flowsheet, p.model, Mocca.stages(p); num_cycles = 1, max_dt = 2.0, first_dt = 0.1)
    kw = (inflow = (:V_feed, :V_lpp), outflow = (:V_product, :V_vacuum), consecutive = 1, min_cycles = 2)
    io = IOBuffer()
    loose = Mocca.simulate_until_mass_balance(p.model, s0, prm, forces, dt; kw..., tol = 1.0, max_cycles = 3, progress = io)
    @test loose.converged && loose.cycles == 2
    lines = split(String(take!(io)), '\n'; keepempty = false)
    @test length(lines) == 3
    @test contains(lines[1], "Mass balance") && contains(lines[1], "1/1 in a row (min 2 cycles)")
    @test contains(lines[end], "converged  below for 1 cycle in a row")
    tight = Mocca.simulate_until_mass_balance(p.model, s0, prm, forces, dt; kw..., tol = 1e-12, max_cycles = 2)
    @test !tight.converged && tight.cycles == 2 && length(tight.history) == 2
    @test tight.history[2] < tight.history[1]

    p = dual_adsorbent_vsa_process(ncells = 10, timestepping = ts)
    @test Mocca.cycle_time(p) ≈ 20.0 + 46.20 + 56.30 + 101.20
    states, steps, s0, prm = run_one_cycle(p)
    check_conservation(p.flowsheet, states, steps, s0, prm, (:SilicaGel, :Zeolite), (:V_feed, :V_lpp), (:V_waste, :V_product, :V_vacuum))
    # The silica gel bed holds back the water
    n_feed = Mocca.stream_totals(states, steps, :V_feed)
    n_link = Mocca.stream_totals(states, steps, :V_link)
    @test n_link[3] < 1e-3 * n_feed[3]
end

@testset "Two-bed VSA with pressure equalisation" begin
    P_H, P_L = 1e5, 0.1e5
    # Both beds start where the cycle has them after equalising, midway between
    # the high and low pressures
    P_eq = (P_H + P_L) / 2
    p = two_bed_vsa_process(ncells = 20, timestepping = TimeStepping(max_dt = 1.0, first_dt = 1.0),
        initial_state = Dict(:A => (pressure = P_eq,), :B => (pressure = P_eq,)))
    @test length(Mocca.stages(p)) == 6
    state0, prm = Mocca.initial_state(p), p.parameters
    r = simulate(p; state0 = state0)
    states, timesteps = r.states, r.timesteps
    @test sum(timesteps) ≈ 80.0

    P_mean(s, b) = sum(s[b][:Pressure]) / length(s[b][:Pressure])
    before = _state_at(states, timesteps, 30.0)
    after = _state_at(states, timesteps, 40.0)
    # A gives gas to B through the equalisation valve, closing most of the gap
    @test after[:V_eq][:MolarFlow][1] > 0
    gap_before = P_mean(before, :A) - P_mean(before, :B)
    gap_after = P_mean(after, :A) - P_mean(after, :B)
    @test gap_before > 0.5 * P_H
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
