using Test
using Mocca
using Jutul

function gradient_test_case(; v_feed = Mocca.HaghpanahConstants{Float64}().v_feed, ncells = 8, kwarg...)
    constants = Mocca.HaghpanahConstants{Float64}(v_feed = v_feed)
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = ncells, kwarg...)
    mm = Mocca.setup_flowsheet_model(fs)
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.setup_process_state(fs[:Bed]; init...))
    prm = Mocca.setup_flowsheet_parameters(fs)
    A = π * constants.r_in^2
    return (; constants, fs, stages, mm, state0, prm, A)
end

# The forces with `change(force_set)` applied to the force set of stage `stage`
function change_stage(forces, stage, change)
    sets = unique(objectid, forces)
    new = IdDict{Any, Any}(f => (i == stage ? change(f) : f) for (i, f) in enumerate(sets))
    return [new[f] for f in forces]
end

@testset "Device force vectorisation" begin
    c = gradient_test_case()
    forces, dt = Mocca.setup_schedule(c.fs, c.mm, c.stages; num_cycles = 1, max_dt = 5.0)
    for f in unique(objectid, forces)
        X, cfg = Jutul.vectorize_forces(f, c.mm)
        g = Jutul.devectorize_forces(f, c.mm, X, cfg)
        @test Jutul.vectorize_forces(g, c.mm)[1] == X
        for d in (:V_feed, :V_product, :V_vacuum)
            @test g[d].law == f[d].law
        end
    end
    # Evacuation: the vacuum valve's conductance and opening, and the ramp
    # (start, stop, λ), temperature and composition of its outlet
    evac = forces[end]
    @test Jutul.vectorization_length(evac[:V_vacuum].outlet, c.mm[:V_vacuum], :outlet, :all) == 3 + 1 + 2
    @test Jutul.vectorization_length(evac[:V_feed].law, c.mm[:V_feed], :law, :all) == 0
end

@testset "Force gradients" begin
    c = gradient_test_case(stage_durations = [5.0, 5.0, 5.0, 5.0])
    forces, dt = Mocca.setup_schedule(c.fs, c.mm, c.stages; num_cycles = 1, max_dt = 1.0)
    G = Mocca.StreamObjective(:V_vacuum; scale = 1e3)
    g = Mocca.force_gradients(c.mm, c.state0, c.prm, forces, dt, G)
    @test length(g.forces) == 4
    function objective(f)
        states, ts, sf = Mocca._simulate_steps(c.mm, c.state0, c.prm, f, dt)
        return Mocca._evaluate_objective(G, c.mm, states, ts, sf)
    end
    @test g.objective ≈ objective(forces)

    # Feed rate in adsorption, a device law
    rate = forces[findfirst(f -> f[:V_feed].law isa Mocca.VolumetricFlow, forces)][:V_feed].law.rate
    function with_rate(q)
        return change_stage(forces, 2, f -> merge(f, Dict(:V_feed => (law = Mocca.VolumetricFlow(q), inlet = f[:V_feed].inlet, outlet = f[:V_feed].outlet))))
    end
    h = 1e-4 * rate
    @test isapprox(g.forces[2][:V_feed].law.rate, (objective(with_rate(rate + h)) - objective(with_rate(rate - h))) / (2h), rtol = 1e-6)

    # Final pressure of the evacuation ramp, a boundary condition
    function with_stop(p)
        return change_stage(forces, 4, function (f)
            d = f[:V_vacuum]
            cond = d.outlet.condition
            r = cond.pressure
            ramp = Mocca.ExponentialRamp(r.start, p, r.λ; t0 = r.t0)
            merge(f, Dict(:V_vacuum => (law = d.law, inlet = d.inlet, outlet = Mocca.PortSetpoint{2}(Mocca.PortCondition(ramp, cond.temperature, cond.composition)))))
        end)
    end
    p = forces[end][:V_vacuum].outlet.condition.pressure.stop
    h = 1e-4 * p
    @test isapprox(g.forces[4][:V_vacuum].outlet.condition.pressure.stop, (objective(with_stop(p + h)) - objective(with_stop(p - h))) / (2h), rtol = 1e-6)
    # A closed device does not depend on its law
    @test g.forces[4][:V_feed].law == Mocca.Closed()

    # Stage durations, with each stage's steps stretched in proportion
    reference = [5.0, 5.0, 5.0, 5.0]
    f_ref, dt_ref = Mocca.setup_schedule(c.fs, c.mm, c.stages; num_cycles = 1, max_dt = 1.0, reference_durations = reference)
    @test dt_ref == dt
    function stretched(i, factor)
        d = copy(reference)
        d[i] *= factor
        s = gradient_test_case(stage_durations = d)
        return Mocca.setup_schedule(s.fs, s.mm, s.stages; num_cycles = 1, max_dt = 1.0, reference_durations = reference)
    end
    function objective_with(f, steps)
        states, ts, sf = Mocca._simulate_steps(c.mm, c.state0, c.prm, f, steps)
        return Mocca._evaluate_objective(G, c.mm, states, ts, sf)
    end
    @test length(g.durations) == 4
    for i in 1:4
        h = 1e-4
        fd = (objective_with(stretched(i, 1 + h)...) - objective_with(stretched(i, 1 - h)...)) / (2h * reference[i])
        @test isapprox(g.durations[i], fd, rtol = 1e-6)
    end
end

@testset "Cyclic steady state gradient" begin
    # Short first steps keep the cycle map smooth through the abrupt start of
    # each stage
    c = gradient_test_case(ncells = 6)
    forces, dt = Mocca.setup_schedule(c.fs, c.mm, c.stages; num_cycles = 1, max_dt = 2.5, first_dt = 0.05)
    loose = Mocca.simulate_to_cyclic_steady_state(c.mm, c.state0, c.prm, forces, dt; tol = 1e-4, max_cycles = 200)
    @test loose.converged
    css = Mocca.newton_cyclic_steady_state(c.mm, loose.state0, c.prm, forces, dt; tol = 1e-11)
    @test css.converged
    @test css.cycles <= 4
    # Each iteration gains at least an order of magnitude
    @test all(css.history[2:end] .< 0.1 .* css.history[1:end-1])

    G = (co2 = Mocca.StreamObjective(:V_vacuum; scale = 1e3), energy = Mocca.VacuumPumpObjective(:V_vacuum; scale = 1e-3))
    g = Mocca.cyclic_steady_state_gradient(c.mm, css.state0, c.prm, forces, dt, G; with_forces = true)
    @test all(r -> r.converged, g)

    # Finite differences, each from a steady state converged by Newton
    reference = [st.duration for st in c.stages]
    function objectives(v_feed, prm; durations = reference)
        s = gradient_test_case(ncells = 6, v_feed = v_feed, stage_durations = durations)
        f, ts = Mocca.setup_schedule(s.fs, s.mm, s.stages; num_cycles = 1, max_dt = 2.5, first_dt = 0.05,
            reference_durations = reference)
        r = Mocca.newton_cyclic_steady_state(s.mm, css.state0, prm, f, ts; tol = 1e-11)
        @test r.converged
        states, steps, sf = Mocca._simulate_steps(s.mm, r.state0, prm, f, ts)
        return map(Gi -> Mocca._evaluate_objective(Gi, s.mm, states, steps, sf), G)
    end
    @test g.co2.objective ≈ objectives(c.constants.v_feed, c.prm).co2

    v = c.constants.v_feed
    h = 1e-4 * v
    jp, jm = objectives(v + h, c.prm), objectives(v - h, c.prm)
    for k in keys(G)
        @test isapprox(g[k].forces[2][:V_feed].law.rate * c.A, (jp[k] - jm[k]) / (2h), rtol = 1e-6)
    end

    hp = 1e-4 .* c.prm[:Bed][:SolidVolume]
    pp = deepcopy(c.prm); pp[:Bed][:SolidVolume] .+= hp
    pm = deepcopy(c.prm); pm[:Bed][:SolidVolume] .-= hp
    jp, jm = objectives(v, pp), objectives(v, pm)
    for k in keys(G)
        @test isapprox(sum(g[k].parameters[:Bed][:SolidVolume] .* hp), (jp[k] - jm[k]) / 2, rtol = 1e-6)
    end

    # Adsorption time
    h = 1e-4 * reference[2]
    jp = objectives(v, c.prm; durations = reference .+ [0, h, 0, 0])
    jm = objectives(v, c.prm; durations = reference .- [0, h, 0, 0])
    for k in keys(G)
        @test isapprox(g[k].durations[2], (jp[k] - jm[k]) / (2h), rtol = 1e-6)
    end
end
