using Test
using Mocca
using Jutul
import Jutul.DictOptimization: DictParameters, free_optimization_parameter!, parameters_gradient

@testset "Cycle metrics" begin
    constants = haghpanah_constants()
    proc = haghpanah_vsa(constants; ncells = 20)
    fs, prm = proc.flowsheet, proc.parameters
    result = simulate(proc; cycles = 2)
    states, ts = result.states, result.timesteps

    w2 = Mocca.cycle_window(ts, 100.0, 2)
    @test sum(ts[w2]) ≈ 100.0
    @test sum(ts[Mocca.cycle_window(ts, 100.0, 1)]) ≈ 100.0

    n_feed = Mocca.stream_totals(states, ts, :V_feed; window = w2)
    n_raff = Mocca.stream_totals(states, ts, :V_product; window = w2)
    n_ext = Mocca.stream_totals(states, ts, :V_vacuum; window = w2)
    # The feed stream has the feed composition
    @test n_feed[1] / sum(n_feed) ≈ constants.y_feed[1] rtol = 1e-6
    # Component balance over the cycle: feed − products = change in the bed
    inv_start = Mocca.component_inventory(fs[:Bed], states[first(w2) - 1][:Bed], prm[:Bed])
    inv_end = Mocca.component_inventory(fs[:Bed], states[last(w2)][:Bed], prm[:Bed])
    @test isapprox(n_feed .- n_raff .- n_ext, inv_end .- inv_start, rtol = 1e-6, atol = 1e-8 * sum(n_feed))

    p = Mocca.purity(n_ext)
    r = Mocca.recovery(n_ext, n_feed)
    @test 0 < p < 1
    @test 0 < r < 1
    @test p > constants.y_feed[1]  # the extract is enriched in CO2
    m_s = Mocca.sorbent_mass(fs[:Bed], prm[:Bed])
    @test m_s ≈ sum(prm[:Bed][:SolidVolume]) * constants.ρ_s
    @test Mocca.productivity(n_ext, 100.0, m_s) > 0

    W = Mocca.vacuum_pump_energy(states, ts, :V_vacuum; window = w2)
    @test W > 0
    e = Mocca.specific_energy_kwh_per_tonne(W, n_ext, constants.molecularMassOfCO2)
    @test isfinite(e) && e > 0
    # The feed valve never runs a vacuum pump
    @test Mocca.vacuum_pump_energy(states, ts, :V_feed; window = w2, discharge_pressure = 1.0) == 0
end

@testset "Cyclic steady state" begin
    p = haghpanah_vsa(; ncells = 10, max_dt = 5.0)
    model, prm = p.model, p.parameters
    state0 = Mocca.initial_state(p)
    forces, dt = Mocca.cycle_schedule(p)

    res = Mocca.simulate_to_cyclic_steady_state(model, state0, prm, forces, dt; tol = 1e-3, max_cycles = 200)
    @test res.converged
    @test res.history[end] < 1e-3
    @test sum(res.timesteps) ≈ 100.0
    @test Mocca.cycle_change(model, res.state0, Mocca.restart_state(model, res.states[end])) ≈ res.history[end]

    # Successive substitution is much slower to get there
    plain = Mocca.simulate_to_cyclic_steady_state(model, state0, prm, forces, dt; tol = 1e-3, max_cycles = res.cycles, acceleration = :none)
    @test !plain.converged
    @test plain.history[end] > res.history[end]

    # At steady state the cycle's feed and products balance
    n_feed = Mocca.stream_totals(res.states, res.timesteps, :V_feed)
    n_out = Mocca.stream_totals(res.states, res.timesteps, :V_product) .+ Mocca.stream_totals(res.states, res.timesteps, :V_vacuum)
    @test isapprox(sum(n_feed), sum(n_out), rtol = 1e-2)
end

@testset "Scale of trace components and empty loadings" begin
    bed = haghpanah_vsa(; ncells = 4).flowsheet[:Bed]
    floors = Mocca.LOADING_SCALE_FRACTION .* Mocca.saturation_capacity(bed.system.isotherm)
    state(y, q) = Mocca.setup_process_state(bed; Pressure = 1e5, Temperature = 298.15, WallTemperature = 298.15,
        y = y, AdsorbedConcentration = q)
    # A bed with nothing adsorbed, at the minimum loading, and the same bed after
    # a change of round-off size: the change is measured against each
    # component's floor
    s0 = state([0.5, 0.5], fill(1e-10, 2, 4))
    s1 = deepcopy(s0)
    s1[:AdsorbedConcentration] .+= 1e-12
    @test Mocca.cycle_change(bed, s0, s1) ≈ 1e-12 / minimum(floors)
    @test Mocca._variable_scale(bed, :AdsorbedConcentration, s0[:AdsorbedConcentration]) == floors
    @test minimum(Mocca._stack_scale(bed, s0)) >= min(maximum(floors), 0.5)
    # A loaded component keeps its own scale, an empty one gets its floor
    q = [fill(2floors[1], 1, 4); fill(1e-10, 1, 4)]
    @test Mocca._variable_scale(bed, :AdsorbedConcentration, q) == [2floors[1], floors[2]]

    # A trace in the gas is measured against its own level, not against the
    # major component
    s0 = state([1e-5, 1 - 1e-5], fill(1.0, 2, 4))
    s1 = deepcopy(s0)
    s1[:y][1, :] .+= 1e-9
    @test Mocca.cycle_change(bed, s0, s1) ≈ 1e-9 / 1e-5
    # and a fraction below the floor against the floor
    s0 = state([1e-10, 1 - 1e-10], fill(1.0, 2, 4))
    s1 = deepcopy(s0)
    s1[:y][1, :] .+= 1e-12
    @test Mocca.cycle_change(bed, s0, s1) ≈ 1e-12 / Mocca.MOLE_FRACTION_SCALE_FLOOR
end

@testset "Flowsheet adjoint gradients" begin
    # Gradients of a product stream with respect to bed parameters, through the
    # devices and cross terms, match central finite differences
    p = haghpanah_vsa(; ncells = 8, durations = (5.0, 5.0, 5.0, 5.0))
    mm, prm = p.model, p.parameters
    state0 = Mocca.initial_state(p)
    forces, dt = Mocca.cycle_schedule(p)
    G(model, state, dt, step_info, forces) = dt * state[:V_vacuum][:MolarFlow][1] * state[:V_vacuum][:InletComposition][1, 1]
    function run(p)
        case = Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = p)
        # Sub-steps must be stored: the adjoint linearises each step that was solved
        sim, cfg = Mocca.setup_process_simulator(mm, state0, p; info_level = -1, nonlinear_tolerance = 1e-8, output_substates = true)
        result = Jutul.simulate!(sim, dt; forces = forces, config = cfg)
        states, ts = Jutul.expand_to_ministeps(result)
        return (case, result, sum(G(mm, s, h, nothing, nothing) for (s, h) in zip(states, ts)))
    end
    case, result, obj = run(prm)
    @test obj > 0
    grad = Jutul.solve_adjoint_sensitivities(case, result, G; info_level = -1)
    # SolidVolume and FluidVolume multiply time differences of the loading and
    # pressure; they need those differences to go through secondary variables
    for p in (:FluidViscosity, :InnerHeatTransferCoeff, :Transmissibilities, :SolidVolume, :FluidVolume)
        h = 1e-4 .* abs.(prm[:Bed][p])
        pp = deepcopy(prm); pp[:Bed][p] .+= h
        pm = deepcopy(prm); pm[:Bed][p] .-= h
        dG_fd = (run(pp)[3] - run(pm)[3]) / 2
        dG_adj = sum(grad[:Bed][p] .* h)
        @test isapprox(dG_adj, dG_fd, rtol = 1e-5)
    end
end

@testset "Column length as an optimisation variable" begin
    # Changing the column length scales the cell lengths, so the gradient runs
    # through the volumes, areas and transmissibilities of the cells, including
    # SolidVolume and FluidVolume. The valves keep a fixed conductance so that
    # only parameters depend on the length, as deps = :parameters requires.
    p = haghpanah_vsa(; ncells = 8, durations = (5.0, 5.0, 5.0, 5.0))
    mm, base, state0 = p.model, p.parameters, Mocca.initial_state(p)
    forces, dt = Mocca.cycle_schedule(p)
    with_length(L) = ((:SolidVolume, :FluidVolume, :CellDx, :WallAreaIn, :WallAreaOut) => L,
        (:Transmissibilities, :DispersionTransmissibilities, :ThermalConductivities) => 1 / L)
    function setup_case(prm, step_info = missing)
        L = prm["L"]
        bed = Dict{Symbol, Any}(k => v for (k, v) in base[:Bed])
        for (names, factor) in with_length(L), k in names
            bed[k] = base[:Bed][k] .* factor
        end
        parameters = merge(base, Dict(:Bed => bed))
        return Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = parameters)
    end
    G(model, state, dt, step_info, forces) = dt * state[:V_vacuum][:MolarFlow][1] * state[:V_vacuum][:InletComposition][1, 1]
    function objective_at(L)
        case = setup_case(Dict("L" => L))
        sim, cfg = Mocca.setup_process_simulator(case.model, case.state0, case.parameters; info_level = -1, nonlinear_tolerance = 1e-8, output_substates = true)
        states, ts = Jutul.expand_to_ministeps(Jutul.simulate!(sim, case.dt; forces = case.forces, config = cfg))
        return sum(G(case.model, s, h, nothing, nothing) for (s, h) in zip(states, ts))
    end
    h = 1e-4
    dG_fd = (objective_at(1.0 + h) - objective_at(1.0 - h)) / (2h)

    case0 = setup_case(Dict("L" => 1.0))
    sim, cfg = Mocca.setup_process_simulator(case0.model, case0.state0, case0.parameters; info_level = -1, nonlinear_tolerance = 1e-8, output_substates = true)
    dprm = DictParameters(Dict("L" => 1.0))
    free_optimization_parameter!(dprm, "L"; abs_min = 0.5, abs_max = 2.0)
    backend_arg = (use_sparsity = true, di_sparse = true, single_step_sparsity = true, do_prep = true, deps = :parameters, deps_ad = :jutul)
    g = parameters_gradient(dprm, G, setup_case; simulator = sim, config = cfg, deps = :parameters, backend_arg = backend_arg)
    @test isapprox(g["L"], dG_fd, rtol = 1e-6)
end

