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
    cell(T, y) = (Temperature = fill(T, 3), y = repeat(y, 1, 3))
    function term(ct, eq, target, dev)
        out = zeros(eq isa typeof(eq_m) ? 2 : 1)
        Jutul.update_cross_term_in_entity!(out, 1, target, target, dev, dev, bed, nothing, ct, eq, 1.0)
        return out
    end
    cp(y) = sum(y .* sys.heat_capacity_gas) * sum(y .* sys.molecular_masses)
    # A link from the top of bed A (cell 3) to the bottom of bed B (cell 1),
    # whose cells differ in temperature and composition
    yA, TA, yB, TB = [0.3, 0.7], 320.0, [0.05, 0.95], 300.0
    A, B = cell(TA, yA), cell(TB, yB)
    for F in (0.02, -0.03)
        dev = device(F, yA, TA, yB, TB)
        y_up, T_up = F > 0 ? (yA, TA) : (yB, TB)
        # What leaves one bed enters the other, component by component and in
        # enthalpy, which is that of the upstream port
        @test term(Mocca.StreamCT(3, 1), eq_m, A, dev) ≈ -term(Mocca.StreamCT(1, 2), eq_m, B, dev)
        @test term(Mocca.StreamCT(3, 1), eq_e, A, dev) ≈ -term(Mocca.StreamCT(1, 2), eq_e, B, dev)
        @test term(Mocca.StreamCT(1, 2), eq_e, B, dev)[1] ≈ -F * cp(y_up) * T_up
    end
    # A stream from a boundary brings its own enthalpy, whatever its pressure
    # and the state of the cell it enters
    feed = [0.15, 0.85]
    e1 = term(Mocca.StreamCT(1, 2), eq_e, cell(300.0, [1e-10, 1 - 1e-10]), device(0.1, feed, 298.15, feed, 298.15))
    e2 = term(Mocca.StreamCT(1, 2), eq_e, cell(350.0, [0.5, 0.5]), device(0.1, feed, 298.15, feed, 298.15; P_in = 0.1e5))
    @test e1 ≈ e2 ≈ [-0.1 * cp(feed) * 298.15]
end

@testset "A port of several cells" begin
    bed = haghpanah_channel_bed(; ncells = 3, channels = 2, thermal = Mocca.WithWall())
    sys = bed.system
    @test Mocca.ports(bed) == (bottom = [1, 4], top = [3, 6])
    @test Mocca.axial_layers(bed) == [[1, 4], [2, 5], [3, 6]]
    @test Mocca.port_weights(bed, [1, 4]) ≈ [0.5, 0.5]
    # Each channel has half the cross-section, so the port conducts as the
    # end of the whole bed
    @test Mocca.port_conductance(bed, :bottom) ≈ Mocca.port_conductance(haghpanah_bed(ncells = 3), :bottom)

    eq_m, eq_e = bed.equations[:mass_conservation], bed.equations[:energy_column]
    dev_model = Mocca.setup_device_model(sys.component_names)
    cp(y) = sum(y .* sys.heat_capacity_gas) * sum(y .* sys.molecular_masses)
    # The top cells differ in temperature and composition
    yt, Tt = ([0.4, 0.6], [0.1, 0.9]), (330.0, 300.0)
    s = (Pressure = [1e5, 1e5, 1.1e5, 1e5, 1e5, 0.9e5], Temperature = [300.0, 300.0, Tt[1], 300.0, 300.0, Tt[2]],
        y = hcat([0.5, 0.5], [0.5, 0.5], yt[1], [0.5, 0.5], [0.5, 0.5], yt[2]))
    top = Mocca.PortStateCT([3, 6], [0.5, 0.5])
    # The device sees the gas of the two cells mixed: the mean pressure and
    # composition, and the temperature that carries their enthalpy
    y_mix = (yt[1] .+ yt[2]) ./ 2
    T_mix = (cp(yt[1]) * Tt[1] + cp(yt[2]) * Tt[2]) / 2 / cp(y_mix)
    F0 = (MolarFlow = [0.0], InletComposition = zeros(2, 1), InletTemperature = [0.0], InletPressure = [0.0],
        OutletComposition = zeros(2, 1), OutletTemperature = [0.0], OutletPressure = [0.0])
    out = zeros(3)
    Jutul.update_cross_term_in_entity!(out, 1, F0, F0, s, s, dev_model, bed, top, nothing, 1.0)
    @test out ≈ -[1e5, T_mix, y_mix[1]]
    @test Mocca.port_state(bed, s, :top).composition ≈ y_mix

    # The device's inlet holds that mixed gas; a stream through it leaves the
    # top cells, each giving up its share of its own gas, and enters another
    # unit at the device's outlet
    device(F) = (MolarFlow = [F], InletComposition = reshape(y_mix, :, 1), InletTemperature = [T_mix], InletPressure = [1e5],
        OutletComposition = reshape([0.15, 0.85], :, 1), OutletTemperature = [298.15], OutletPressure = [0.5e5])
    function terms(ct, eq, dev)
        n = eq === eq_m ? 2 : 1
        return [(o = zeros(n); Jutul.update_cross_term_in_entity!(o, i, s, s, dev, dev, bed, nothing, ct, eq, 1.0); o)
            for i in eachindex(ct.cells)]
    end
    out_of_top = Mocca.StreamCT([3, 6], [0.5, 0.5], 1)
    F = 0.02
    m = terms(out_of_top, eq_m, device(F))
    @test m[1] ≈ 0.5 * F .* yt[1] && m[2] ≈ 0.5 * F .* yt[2]
    @test sum(m) ≈ F .* y_mix
    e = terms(out_of_top, eq_e, device(F))
    @test sum(e)[1] ≈ F * cp(y_mix) * T_mix
    # Flowing back, the top cells share the gas of the device's outlet
    m = terms(out_of_top, eq_m, device(-F))
    @test m[1] ≈ m[2] ≈ -0.5 * F .* [0.15, 0.85]
end

@testset "A bed in equal channels is the bed" begin
    # Split side by side into two equal channels, a bed fed evenly across its
    # end behaves as one channel
    one = haghpanah_vsa(; ncells = 10, thermal = Mocca.Adiabatic())
    two = haghpanah_vsa(; ncells = 10, thermal = Mocca.Adiabatic(),
        bed_model = haghpanah_channel_bed(; ncells = 10, channels = 2))
    r1, r2 = simulate(one), simulate(two)
    @test length(r1.states) == length(r2.states)
    bed2 = two.flowsheet[:Bed]
    s1, s2 = r1.states[end], r2.states[end]
    for k in (:Pressure, :Temperature, :y, :AdsorbedConcentration)
        @test Mocca.axial_profile(bed2, s2[:Bed][k]) ≈ s1[:Bed][k] rtol = 1e-6
    end
    # Both channels are the same
    @test s2[:Bed][:y][:, 1:10] ≈ s2[:Bed][:y][:, 11:20] rtol = 1e-8
    for d in (:V_feed, :V_product, :V_vacuum)
        @test Mocca.stream_totals(r2.states, r2.timesteps, d) ≈ Mocca.stream_totals(r1.states, r1.timesteps, d) rtol = 1e-6
    end
    @test performance(r2).recovery ≈ performance(r1).recovery rtol = 1e-6
    # Plots along the bed show the cross-section averages
    @test Mocca.plot_state(s2[:Bed], bed2) isa Mocca.Plotting.Figure
    @test Mocca.plot_profiles(r2.states, r2.timesteps, bed2, [30.0, 100.0]; unit = :Bed) isa Mocca.Plotting.Figure
end

@testset "Channels that differ conserve moles" begin
    # The channels hold different amounts and do not exchange gas or heat
    # along the bed, so their gas differs, and when the bed is evacuated
    # through its bottom each bottom cell gives up its own gas
    p = haghpanah_vsa(; ncells = 10, thermal = Mocca.Adiabatic(),
        bed_model = haghpanah_channel_bed(; ncells = 10, channels = 2, porosity = [0.3, 0.45]))
    bed = p.flowsheet[:Bed]
    N = bed.data_domain[:neighbors]
    between = [f for f in axes(N, 2) if abs(N[1, f] - N[2, f]) == 10]
    for k in (:Transmissibilities, :DispersionTransmissibilities, :ThermalConductivities)
        p.parameters[:Bed][k][between] .= 0.0
    end
    state0, prm = Mocca.initial_state(p), p.parameters
    r = simulate(p; state0 = state0)
    @test maximum(abs(s[:Bed][:y][1, 1] - s[:Bed][:y][1, 11]) for s in r.states) > 0.01
    held(s) = Mocca.component_inventory(bed, s[:Bed], prm[:Bed])
    n_in = Mocca.stream_totals(r.states, r.timesteps, :V_feed)
    n_out = Mocca.stream_totals(r.states, r.timesteps, :V_product) .+ Mocca.stream_totals(r.states, r.timesteps, :V_vacuum)
    @test held(r.states[end]) .- held(state0) ≈ n_in .- n_out rtol = 1e-7 atol = 1e-8 * sum(n_in)
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
