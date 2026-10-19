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
    V_s = Mocca.adsorbent_volume(fs[:Bed], prm[:Bed])
    @test V_s ≈ π * constants.r_in^2 * constants.L * (1 - constants.Φ)
    @test Mocca.productivity(n_ext, 100.0, V_s) ≈ n_ext[1] / (V_s * 100.0)

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

    io = IOBuffer()
    res = Mocca.simulate_to_cyclic_steady_state(model, state0, prm, forces, dt; tol = 1e-3, max_cycles = 200, progress = io)
    @test res.converged
    @test res.history[end] < 1e-3
    # One line per cycle, then the end of the phase
    lines = split(String(take!(io)), '\n'; keepempty = false)
    @test length(lines) == res.cycles + 1
    @test all(contains("Anderson"), lines)
    @test contains(lines[1], "cycle   1") && contains(lines[1], "change ")
    @test contains(lines[end], "$(res.cycles) cycles") && contains(lines[end], "< 1.0e-03") && contains(lines[end], "  converged")
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

@testset "Cyclic steady state progress" begin
    O = Mocca.Operation
    @test O._sci(2.346e-3) == "2.35e-03"
    @test O._sci(1e-9) == "1.0e-09"
    @test O._sci(9.996e-4) == "1.0e-03"
    @test O._sci(0.0) == "0"
    @test O._duration(4.12) == "4.1 s"
    @test O._duration(32.4) == "32 s"
    @test O._duration(75.0) == "1 min 15 s"
    @test O._duration(3700.0) == "1 h 2 min"
    @test O._bar(0.5; width = 4) == "▕██  ▏"
    @test O._bar(0.5625; width = 4) == "▕██▎ ▏"
    @test O._bar(2.0; width = 4) == "▕████▏"
    # Two of four decades from the first change to the tolerance, by the
    # smallest change so far
    @test O._log_fraction([1e-1, 1e-3], 1e-5) ≈ 0.5
    @test O._log_fraction([1e-1, 1e-3, 1e-2], 1e-5) ≈ 0.5
    @test O._log_fraction([1e-1, 1.0], 1e-5) == 0
    # At a tenfold drop per cycle, 3e-5 is reached in two more cycles; the
    # first cycle, from the initial guess, does not count
    @test O._cycles_to_reach([1e3, 1e-1, 1e-2, 1e-3], 3e-5) == 2
    @test O._cycles_to_reach([1e3, 1e-1, 1e-2, 1e-1, 1e-3], 3e-5) == 3
    @test isnothing(O._cycles_to_reach([1e3, 1e-1, 1e-2], 3e-5))
    @test isnothing(O._cycles_to_reach([1e3, 1e-2, 2e-2, 1e-2], 1e-4))
    # Four more cycles below the tolerance, unless the minimum is further away
    @test O._mass_balance_cycles_left([0.1, 0.01, 0.001], 0.005, 5, 0) == (4, 1)
    @test O._mass_balance_cycles_left([0.1, 0.01, 0.001], 0.005, 5, 50) == (47, 1)
    @test O._mass_balance_cycles_left([1.0, 0.1, 0.03, 0.009], 0.005, 5, 0) == (5, 0)
    @test O._mass_balance_cycles_left([0.1, 0.03], 0.005, 5, 10) == (8, 0)

    # What a terminal shows of a live printout, which only returns to the start
    # of a line, starts a new one, clears to its end and moves up
    function terminal_screen(out)
        rows, r, col = [""], 1, 0
        for m in eachmatch(r"\r|\n|\e\[K|\e\[(\d+)A|[^\r\n\e]+", out)
            t = m.match
            if t == "\r"
                col = 0
            elseif t == "\n"
                r += 1
                r > length(rows) && push!(rows, "")
                col = 0
            elseif t == "\e[K"
                rows[r] = first(rows[r], col)
            elseif startswith(t, "\e[")
                r -= parse(Int, m.captures[1])
            else
                rows[r] = first(rows[r], col) * t
                col += length(t)
            end
        end
        return rows
    end

    # A live printout of a single phase rewrites one line in place, kept within
    # the display
    io = IOBuffer()
    pr = Mocca.CycleProgress(IOContext(io, :displaysize => (24, 60)); live = true)
    O.start_phase!(pr, "Anderson"; tol = 1e-4)
    O.report_cycle!(pr, 1, 1e-1)
    O.report_cycle!(pr, 2, 1e-2)
    O.finish_phase!(pr, 2, false)
    O.finish_progress!(pr, false)
    screen = terminal_screen(String(take!(io)))
    @test length(screen) == 3 && screen[end] == ""
    @test all(l -> textwidth(l) <= 59, screen)
    @test occursin(r"^  Anderson +2 cycles +change 1\.0e-02 ≥ 1\.0e-04 ", screen[1])
    @test startswith(screen[2], "Not converged after 2 cycles simulated")

    # With two phases planned, each has its own line and bar: the second
    # waits for the first, which is done rather than converged when it hands
    # over
    io = IOBuffer()
    out = ""
    show_screen() = (out *= String(take!(io)); terminal_screen(out))
    pr = Mocca.CycleProgress(IOContext(io, :displaysize => (24, 200)); live = true)
    O.plan_phases!(pr, [("Anderson", 1e-3, "cycle"), ("Newton", 1e-9, "iteration")])
    O.start_phase!(pr, "Anderson"; tol = 1e-3)
    O.report_cycle!(pr, 1, 1e-1)
    screen = show_screen()
    @test length(screen) == 2
    @test contains(screen[1], "Anderson") && contains(screen[1], "cycle   1") && contains(screen[1], "change 1.0e-01 → 1.0e-03")
    @test contains(screen[2], "Newton") && contains(screen[2], "to 1.0e-09") && contains(screen[2], "  0%  waiting for Anderson")
    # The bars line up
    @test findfirst(==('▕'), collect(screen[1])) == findfirst(==('▕'), collect(screen[2]))
    O.report_cycle!(pr, 2, 1e-4)
    O.finish_phase!(pr, 2, true)
    O.start_phase!(pr, "Newton"; tol = 1e-9, noun = "iteration", estimate = false)
    O.report_cycle!(pr, 0, 1e-4)
    O.report_cycle!(pr, 1, 1e-10)
    O.finish_phase!(pr, 1, true)
    O.finish_progress!(pr, true)
    screen = show_screen()
    @test length(screen) == 4 && screen[end] == ""
    @test contains(screen[1], "2 cycles") && contains(screen[1], "100%  done, on to Newton") && !contains(screen[1], "converged")
    @test contains(screen[2], "1 iteration ") && contains(screen[2], "< 1.0e-09") && contains(screen[2], "100%  converged")
    @test startswith(screen[3], "Converged after 4 cycles simulated")

    # Without a terminal, each cycle gets its own line and the waiting phase
    # none
    io = IOBuffer()
    pr = Mocca.CycleProgress(io; live = false)
    O.plan_phases!(pr, [("Anderson", 1e-3, "cycle"), ("Newton", 1e-9, "iteration")])
    O.start_phase!(pr, "Anderson"; tol = 1e-3)
    O.report_cycle!(pr, 1, 1e-1)
    O.report_cycle!(pr, 2, 1e-4)
    O.finish_phase!(pr, 2, true)
    out = String(take!(io))
    lines = split(out, '\n'; keepempty = false)
    @test length(lines) == 3 && all(contains("Anderson"), lines)
    @test !contains(out, "\e") && !any(startswith("  Newton"), lines)
    @test contains(lines[end], "done, on to Newton")

    # A phase without a target, such as a set number of plain cycles, shows
    # its change and a bar of the cycles done
    io = IOBuffer()
    pr = Mocca.CycleProgress(io; live = false)
    O.plan_phases!(pr, [("Cycles", NaN, "cycle"), ("Newton", 1e-9, "iteration")])
    O.start_phase!(pr, "Cycles"; tol = NaN, estimate = false)
    O.report_cycle!(pr, 1, 0.3; fraction = 1 / 3)
    O.finish_phase!(pr, 3, true)
    lines = split(String(take!(io)), '\n'; keepempty = false)
    @test contains(lines[1], "change 3.0e-01 ") && contains(lines[1], " 33%") && !contains(lines[1], "→")
    @test contains(lines[end], "3 cycles") && contains(lines[end], "done, on to Newton")
    @test !any(contains("NaN"), lines)
    @test isnothing(O.cycle_progress(false))
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
    @test Mocca.Operation._variable_scale(bed, :AdsorbedConcentration, s0[:AdsorbedConcentration]) == floors
    @test minimum(Mocca.Operation._stack_scale(bed, s0)) >= min(maximum(floors), 0.5)
    # A loaded component keeps its own scale, an empty one gets its floor
    q = [fill(2floors[1], 1, 4); fill(1e-10, 1, 4)]
    @test Mocca.Operation._variable_scale(bed, :AdsorbedConcentration, q) == [2floors[1], floors[2]]

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

