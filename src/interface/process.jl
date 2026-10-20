# The Process: one object holding a flowsheet, its cycle, the initial state and
# solver options, with simulation, cyclic steady state and performance.

"""
Options of a [`Process`](@ref) that set how its performance is measured.
"""
struct ProcessOptions
    key_component::Int
    ambient_pressure::Float64
    pump_efficiency::Float64
    heat_capacity_ratio::Float64
    ramp_rate::Float64
end

"Last cyclic steady state and cycle Jacobian, to start the next solve from."
mutable struct WarmStart
    state::Any
    jacobian::Any
end
WarmStart() = WarmStart(nothing, nothing)

"""
    Process(fs, cycle; name = "", timestepping = TimeStepping(), key_component = 1,
        ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = 1.4,
        ramp_rate = 0.5, initial_state = Dict(), parameters = Dict())

A cyclic adsorption process: the flowsheet `fs` (built with
[`Flowsheet(gas, feed)`](@ref Flowsheet), [`add_bed!`](@ref) and
[`add_device!`](@ref)) operated by `cycle` (see [`Cycle`](@ref)).

- `timestepping`: see [`TimeStepping`](@ref).
- `key_component`: the component to capture, by index or name, for
  [`performance`](@ref).
- `ambient_pressure` [Pa]: vacuum pumps discharge to it and feed compressors
  take gas from it; `pump_efficiency` and `heat_capacity_ratio` are those of
  both.
- `ramp_rate` [1/s]: of the pressure ramps of steps that do not set their own.
- `initial_state`: per unit with holdup, a NamedTuple with any of
  `pressure`, `temperature`, `wall_temperature` and `composition`. By default
  each unit starts at the lowest pressure of its cycle (the feed pressure if
  no step sets it), at the feed temperature, full of light product, with the
  sorbent of a bed at equilibrium (see [`setup_unit_state`](@ref)).
- `parameters`: per unit, values of its Jutul parameters to override, such as
  `Dict(:Bed => Dict(:AdsorbentDensity => 1200.0))`.

Simulate it with [`simulate`](@ref), find its cyclic steady state with
[`cyclic_steady_state`](@ref), and measure it with [`performance`](@ref).

The process keeps its own copies of `fs` and `cycle`, so changing them
afterwards does not change the process, and several processes can be built
from the same ones. Change a process through its design variables (see
[`add_variable!`](@ref)); `p.flowsheet` and `p.cycle` are its copies.
"""
mutable struct Process
    name::String
    flowsheet::Flowsheet
    cycle::Cycle
    timestepping::TimeStepping
    options::ProcessOptions
    initial_conditions::Dict{Symbol, Any}
    parameter_settings::Dict{Symbol, Dict{Symbol, Any}}
    variables::Vector{Any}
    model::Jutul.MultiModel
    parameters::Dict{Symbol, Any}
    compiled::Union{Nothing, CompiledCycle}
    reference_durations::Union{Nothing, Vector{Float64}}
    warm::WarmStart
end

function Process(fs::Flowsheet, cycle::Cycle; name = "", timestepping = TimeStepping(), key_component = 1,
        ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = 1.4, ramp_rate = 0.5,
        initial_state = Dict{Symbol, Any}(), parameters = Dict{Symbol, Any}())
    # The process sets up and later changes its flowsheet (dispersion
    # velocities, default boundaries, the feed) and its steps (design
    # variables), so it works on copies of its own
    fs = deepcopy(fs)
    cycle = deepcopy(cycle)
    k = key_component isa Integer ? Int(key_component) : component_index(fs.gas, key_component)
    1 <= k <= number_of_components(fs.gas) || input_error("process", "key_component", "the gas has $(number_of_components(fs.gas)) components, got $k")
    options = ProcessOptions(k, ambient_pressure, pump_efficiency, heat_capacity_ratio, ramp_rate)
    compiled = _compile!(fs, cycle, ramp_rate)
    _set_dispersion_velocities!(fs, compiled)
    for u in holdup_units(fs)
        isempty(vcat([devices_at(fs, u, p) for p in keys(ports(fs[u]))]...)) && @warn "Unit $u has no devices attached"
    end
    model = setup_flowsheet_model(fs)
    ics = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in pairs(initial_state))
    for u in keys(ics)
        u in holdup_units(fs) || input_error("initial_state", String(u), "$u is not a unit with holdup of the flowsheet")
    end
    settings = Dict{Symbol, Dict{Symbol, Any}}(Symbol(u) => Dict{Symbol, Any}(Symbol(n) => v for (n, v) in pairs(d)) for (u, d) in pairs(parameters))
    p = Process(String(name), fs, cycle, timestepping, options, ics, settings, Any[], model, Dict{Symbol, Any}(), compiled, nothing, WarmStart())
    p.parameters = _process_parameters(p)
    return p
end

function _process_parameters(p::Process)
    fs = p.flowsheet
    prm = setup_flowsheet_parameters(fs)
    for b in beds(fs)
        get(fs.specs, b, nothing) isa BedSpec && (prm[b][:AxialDispersion] .= bed_dispersion(fs, b))
    end
    for (u, d) in p.parameter_settings
        haskey(prm, u) || input_error("parameters", String(u), "$u is not a unit of the flowsheet")
        for (n, v) in d
            haskey(prm[u], n) || input_error("parameters of $u", String(n), "unknown parameter; $u has $(join(sort(collect(keys(prm[u]))), ", "))")
            prm[u][n] .= v
        end
    end
    return prm
end

# Beds whose dispersion follows the adsorption velocity take that velocity;
# their AxialDispersion parameter is set from it in _process_parameters
function _set_dispersion_velocities!(fs, compiled)
    all_v = vcat(values(compiled.velocities)...)
    for b in beds(fs)
        spec = get(fs.specs, b, nothing)
        spec isa BedSpec && spec.dispersion_velocity == :adsorption && isnothing(spec.dispersion) || continue
        v = haskey(compiled.velocities, b) ? maximum(compiled.velocities[b]) : (isempty(all_v) ? nothing : maximum(all_v))
        if isnothing(v)
            @warn "No adsorption step sets a feed velocity, so the axial dispersion of bed $b uses 1 m/s; pass dispersion_velocity to add_bed! to set it"
            continue
        end
        fs.specs[b].velocity_used = v
    end
end

# Default boundaries follow the feed and the cycle: they take the feed
# temperature and the composition of their role, and extract and waste lines
# sit at the lowest pressure of their bed, the others at the feed pressure.
# Without `compiled`, the pressures are left as they are.
function _refresh_default_boundaries!(fs, compiled = nothing)
    for r in fs.default_boundaries
        role = device_role(fs, r.name)
        if isnothing(compiled)
            P = boundary_condition(fs, r.name, r.port).pressure
        else
            P = fs.feed.pressure
            if role in (:extract, :waste)
                b = first(attachments(fs, r.name)).unit.name
                haskey(compiled.pressures, b) && (P = minimum(compiled.pressures[b]))
            end
        end
        # On the network, so that the boundary stays a default one
        set_boundary!(fs.network, r.name, r.port, PortCondition(pressure = P, temperature = fs.feed.temperature, composition = role_composition(fs, role)))
    end
end

# Compile the cycle and update the default boundaries. The steps take the
# temperature and composition of a boundary from the network, so those are
# updated first; the pressures follow from the compiled cycle.
function _compile!(fs, cycle, ramp_rate)
    _refresh_default_boundaries!(fs)
    compiled = compile_cycle(cycle, fs; ramp_rate = ramp_rate)
    _refresh_default_boundaries!(fs, compiled)
    return compiled
end

"Recompile the cycle after a change to its steps or the feed."
function _refresh!(p::Process)
    p.compiled = _compile!(p.flowsheet, p.cycle, p.options.ramp_rate)
    return p
end

_compiled(p::Process) = isnothing(p.compiled) ? _refresh!(p).compiled : p.compiled

"""
    stages(p::Process) → Vector{Stage}

The stages of the process's cycle, with the device settings of each.
"""
stages(p::Process) = _compiled(p).stages

"Duration [s] of one cycle."
cycle_time(p::Process) = sum(st -> st.duration, stages(p))

"""
    reset_step_pattern!(p)

Take the current stage durations as the pattern of time steps that changed
durations stretch (see [`TimeStepping`](@ref)).
"""
reset_step_pattern!(p::Process) = (p.reference_durations = nothing; p)

function _keeps_step_pattern(p::Process)
    k = p.timestepping.keep_step_pattern
    k isa Bool && return k
    return any(v -> any(pl -> pl isa StepPlace && pl.field == :duration, v.places), p.variables)
end

function _reference_durations(p::Process)
    _keeps_step_pattern(p) || return nothing
    d = [st.duration for st in stages(p)]
    if isnothing(p.reference_durations) || length(p.reference_durations) != length(d)
        p.reference_durations = d
    end
    return p.reference_durations
end

"""
    cycle_schedule(p; cycles = 1) → (forces, timesteps)

Forces and time steps for `cycles` cycles of the process; see
[`setup_schedule`](@ref).
"""
function cycle_schedule(p::Process; cycles = 1)
    ts = p.timestepping
    return setup_schedule(p.flowsheet, p.model, stages(p); num_cycles = cycles, max_dt = ts.max_dt,
        first_dt = ts.first_dt, reference_durations = _reference_durations(p))
end

"""
    initial_state(p) → Dict

The initial state of every unit; see `initial_state` in [`Process`](@ref).
"""
function initial_state(p::Process)
    fs = p.flowsheet
    compiled = _compiled(p)
    unit_states = Dict{Symbol, Any}()
    for u in holdup_units(fs)
        ic = get(p.initial_conditions, u, NamedTuple())
        P = haskey(ic, :pressure) ? ic.pressure : (haskey(compiled.pressures, u) ? minimum(compiled.pressures[u]) : fs.feed.pressure)
        T = get(ic, :temperature, fs.feed.temperature)
        spec = get(fs.specs, u, nothing)
        Tw = get(ic, :wall_temperature, spec isa BedSpec ? spec.column.ambient_temperature : T)
        y = haskey(ic, :composition) ? composition(fs.gas, ic.composition; location = "initial state of $u") : fs.feed.light_product
        unit_states[u] = setup_unit_state(fs[u]; pressure = P, temperature = T, composition = y, wall_temperature = Tw)
    end
    return setup_flowsheet_state(fs; unit_states...)
end

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

"""
Result of simulating a [`Process`](@ref): the `states` and `timesteps` of
every step the solver took, the number of `cycles` they span, the initial
state `state0` of the last of them, and for a cyclic steady state whether it
`converged`, the change per iteration or cycle in `history` and the cycle
Jacobian if Newton's method found it.
"""
struct ProcessResult
    process::Process
    states::Vector{Any}
    timesteps::Vector{Float64}
    cycles::Int
    state0::Any
    converged::Union{Nothing, Bool}
    history::Vector{Float64}
    iterations::Int
    jacobian::Any
end

"""
    simulate(p::Process; cycles = 1, state0 = initial_state(p), info_level = -1, kwarg...) → ProcessResult

Simulate `cycles` cycles of the process. Extra keywords go to
`simulate_process`.
"""
function Jutul.simulate(p::Process; cycles = 1, state0 = nothing, info_level = -1, kwarg...)
    forces, dt = cycle_schedule(p; cycles = cycles)
    x0 = isnothing(state0) ? initial_state(p) : state0
    case = MoccaCase(p.model, dt, forces; state0 = x0, parameters = p.parameters)
    states, ts = simulate_process(case; output_substates = true, info_level = info_level, kwarg...)
    isapprox(sum(ts), sum(dt)) || throw(SimulationFailure("the simulation stopped at t = $(sum(ts)) s of $(sum(dt)) s"))
    n = cycles
    if cycles > 1
        w = cycle_window(ts, cycle_time(p), cycles)
        last0 = restart_state(p.model, states[first(w) - 1])
    else
        last0 = x0
    end
    return ProcessResult(p, collect(Any, states), ts, n, last0, nothing, Float64[], 0, nothing)
end

# Newton's method starts after this many plain cycles. Each time it does not
# converge, Anderson acceleration goes on to the next of these cycle changes,
# and Newton's method tries again.
const DEFAULT_HANDOVER = 3
const FALLBACK_HANDOVERS = (1e-1, 1e-3)

# Beds with more cells than this start from the steady state on half as many
const DEFAULT_COARSE_CELLS = 40

"""
    cyclic_steady_state(p; method = :newton, tol = nothing, max_cycles = 500,
        warm_start = true, handover = 3, krylov = :auto, coarse_cells = 40,
        progress = isinteractive(), info_level = -1, kwarg...) → ProcessResult

The cyclic steady state of the process, with the states of its last cycle:

- `:newton` (default): `handover` plain cycles (3), each from the end of the
  last, then Newton's method on the cycle map (see
  [`newton_cyclic_steady_state`](@ref)) to `tol` (default `1e-9`). The plain
  cycles let the fastest modes die out, and Newton's method controls its
  steps, so it needs no closer start. A `handover` below 1 is a cycle change
  instead: Anderson acceleration until the change is below it. If Newton's
  method fails or does not converge, Anderson acceleration goes on to a
  change of `1e-1`, and then if need be `1e-3`, and Newton's method tries
  again from each. With `warm_start`, it starts from the last
  steady state found for this process, and its Jacobian, which is much faster
  when the settings have changed a little. With `krylov = true`, each Newton
  step is solved by GMRES instead of from the full Jacobian, which pays off
  for fine grids, where the variables far outnumber the slow modes of the
  cycle; `false` always builds the full Jacobian. The default, `:auto`,
  decides on each grid: GMRES when the beds have more than 800 variables
  (a bed of about 130 cells with two components, 100 with three), the full
  Jacobian below; `krylov = n` sets that number of variables to `n`.

  On a fine grid, with a bed of more than `coarse_cells` cells, the steady
  state is found first on beds of half as many cells (and so on, down to
  `coarse_cells` or fewer), each interpolated onto the next finer grid as the
  start of Newton's method there, which then needs only a few iterations; a
  coarser grid is converged to `1e-6` only. Starting fine grids cold is slow,
  and on some processes (such as the wet LPP VSA from about 80 cells)
  Newton's method does not converge from a few plain cycles. The coarser
  grids only give the start: the steady state is that of `p`. If Newton's
  method does not converge from the coarser grid's steady state, `p` is
  started from its initial state as above. Pass `coarse_cells = nothing` to
  start from the initial state on any grid.
- `:anderson`: [`simulate_to_cyclic_steady_state`](@ref) to `tol` (default
  `1e-6`).
- `:mass_balance`: [`simulate_until_mass_balance`](@ref) on the key component,
  between the feed and light inlet devices and the outlet devices.

With `progress`, on by default in an interactive session such as the REPL,
each cycle is printed as it finishes: how many have run, how much the state
still changes over a cycle against `tol`, a bar of how far it has come, and an
estimate of the cycles and time left. In a terminal the line is rewritten in
place, so only the end of each phase stays on screen. `progress` can also be
an `IO` to print to, or a [`CycleProgress`](@ref).
"""
function cyclic_steady_state(p::Process; method = :newton, tol = nothing, max_cycles = 500, warm_start = true,
        handover = DEFAULT_HANDOVER, krylov = :auto, coarse_cells = DEFAULT_COARSE_CELLS, progress = isinteractive(),
        info_level = -1, kwarg...)
    forces, dt = cycle_schedule(p)
    model, prm = p.model, p.parameters
    has_warm = warm_start && !isnothing(p.warm.state)
    x0 = has_warm ? restart_state(model, p.warm.state) : initial_state(p)
    pr = cycle_progress(progress; info_level = info_level)
    name = isempty(p.name) ? "the process" : p.name
    if method == :newton
        tol = something(tol, 1e-9)
        (handover isa Integer ? handover >= 0 : 0 < handover < 1) ||
            input_error("cyclic_steady_state", "handover", "give a number of cycles, or a cycle change below 1, got $handover")
        (krylov isa Integer || krylov === :auto) ||
            input_error("cyclic_steady_state", "krylov", "use true, false, :auto or a number of variables, got $(repr(krylov))")
        grids = has_warm ? Process[] : _coarse_grids(p, coarse_cells)
        start = "$(_handover_text(handover)), then Newton's method"
        middle = [_cells_text(q) for q in grids[2:end]]
        report_note!(pr, "Cyclic steady state of $name: " * (has_warm ? "Newton's method to $(_sci(tol)) from the last steady state" :
            isempty(grids) ? "$start to $(_sci(tol))" :
            "$start on $(_cells_text(grids[1])) to $(_sci(max(tol, COARSE_TOL))); from there Newton's method on " *
            (isempty(middle) ? "" : "$(join(middle, ", ", " and ")) to $(_sci(max(tol, COARSE_TOL))) and on ") *
            "$(_cells_text(p)) to $(_sci(tol))"))
        p.timestepping.first_dt < p.timestepping.max_dt ||
            @warn "Newton's method needs a smooth cycle map; start each stage with short steps (TimeStepping(first_dt = 0.05))" maxlog = 1
        s = _newton_steady_state(p, pr; tol = tol, warm = has_warm, handover = handover, krylov = krylov,
            coarse_cells = coarse_cells, max_cycles = max_cycles, info_level = info_level, kwarg...)
        r = s.result
        finish_progress!(pr, r.converged)
        r.converged || @warn "Cyclic steady state not converged by Newton's method" r.history
        # Only a steady state is a good start for the next solve; otherwise the
        # last one found, if any, stays
        r.converged && ((p.warm.state, p.warm.jacobian) = (r.state0, r.jacobian))
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, s.history, s.cycles, r.jacobian)
    elseif method == :anderson
        tol = something(tol, 1e-6)
        report_note!(pr, "Cyclic steady state of $name: Anderson acceleration to $(_sci(tol))" * (has_warm ? " from the last steady state" : ""))
        r = simulate_to_cyclic_steady_state(model, x0, prm, forces, dt; tol = tol, max_cycles = max_cycles,
            progress = pr, info_level = info_level, kwarg...)
        finish_progress!(pr, r.converged)
        r.converged && (p.warm.state = r.state0)
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, r.history, r.cycles, nothing)
    elseif method == :mass_balance
        fs = p.flowsheet
        inflow = Tuple(vcat(devices_with_role(fs, :feed), devices_with_role(fs, :light_inlet)))
        outflow = Tuple(vcat([devices_with_role(fs, r) for r in OUTLET_ROLES]...))
        kw = isnothing(tol) ? NamedTuple() : (tol = tol,)
        k = p.options.key_component
        report_note!(pr, "Cyclic steady state of $name: $(fs.gas.names[k]) mass balance" * (has_warm ? " from the last steady state" : ""))
        r = simulate_until_mass_balance(model, x0, prm, forces, dt; inflow = inflow, outflow = outflow,
            component = k, max_cycles = max_cycles, progress = pr, info_level = info_level, kw..., kwarg...)
        finish_progress!(pr, r.converged)
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, r.history, r.cycles, nothing)
    else
        input_error("cyclic_steady_state", "method", "use :newton, :anderson or :mass_balance, got $method")
    end
end

# Cycle change to which a coarser grid is converged when it only gives the
# start for a finer one: its interpolation onto the finer grid is far rougher
const COARSE_TOL = 1e-6

# Newton's method for the cyclic steady state of `p`, as cyclic_steady_state
# runs it: from the last steady state if `warm`; otherwise, if the beds have
# more than `coarse_cells` cells, from the steady state on half as many,
# interpolated (found the same way, so possibly from coarser grids still);
# otherwise, or if that fails, from the initial state, after `handover`
# plain cycles or Anderson acceleration, with Anderson acceleration to
# tighter hand-overs as the fallback. Returns the result of the last Newton
# solve, the change of every cycle or iteration and the cycles simulated.
function _newton_steady_state(p::Process, pr; tol, warm = false, handover = DEFAULT_HANDOVER, krylov = :auto,
        coarse_cells = DEFAULT_COARSE_CELLS, max_cycles = 500, info_level = -1, in_sequence = false, kwarg...)
    forces, dt = cycle_schedule(p)
    model, prm = p.model, p.parameters
    history = Float64[]
    cycles = 0
    if warm
        plan_phases!(pr, [("Newton", tol, "iteration")])
        r = try
            newton_cyclic_steady_state(model, restart_state(model, p.warm.state), prm, forces, dt; tol = tol, max_iterations = 40,
                jacobian = p.warm.jacobian, krylov = krylov, progress = pr, info_level = info_level, kwarg...)
        catch e
            # A start too far from the new steady state can make the
            # solver fail; anything else is a bug and is passed on
            _is_solver_failure(e) || rethrow()
            report_note!(pr, "  Newton's method from the last steady state failed ($(sprint(showerror, e)))")
            @debug "Warm-started Newton's method failed" exception = (e, catch_backtrace())
            nothing
        end
        !isnothing(r) && r.converged && return (result = r, history = r.history, cycles = r.cycles)
        if !isnothing(r)
            append!(history, r.history)
            cycles += r.cycles
        end
        report_note!(pr, "  Starting again from the initial state")
    end

    # From a coarser grid
    q = _coarse_copy(p, coarse_cells)
    if !isnothing(q)
        c = _newton_steady_state(q, pr; tol = max(tol, COARSE_TOL), handover = handover, krylov = krylov,
            coarse_cells = coarse_cells, max_cycles = max_cycles, info_level = info_level, in_sequence = true, kwarg...)
        append!(history, c.history)
        cycles += c.cycles
        if c.result.converged
            report_note!(pr, "  $(_cells_text(p)): Newton's method from the steady state on $(_cells_text(q))")
            plan_phases!(pr, [("Newton", tol, "iteration")])
            r = try
                newton_cyclic_steady_state(model, _interpolate_state(p, c.result.state0), prm, forces, dt; tol = tol,
                    max_iterations = 40, krylov = krylov, progress = pr, info_level = info_level, kwarg...)
            catch e
                _is_solver_failure(e) || rethrow()
                @debug "Newton's method from the coarser grid failed" exception = (e, catch_backtrace())
                nothing
            end
            if !isnothing(r)
                append!(history, r.history)
                cycles += r.cycles
                r.converged && return (result = r, history = history, cycles = cycles)
            end
        end
        report_note!(pr, "  Newton's method did not converge from the coarser grid; $(_cells_text(p)) from the initial state")
    end

    # From the initial state: plain cycles or Anderson acceleration to
    # `handover`, then Newton's method. If Newton's method fails or does not
    # converge from there, Anderson acceleration goes on, from whichever start
    # changed least, to the next tighter hand-over, and Newton's method tries
    # again.
    in_sequence && isnothing(q) && report_note!(pr, "  $(_cells_text(p)): $(_handover_text(handover)), then Newton's method")
    stages = Any[handover; [h for h in FALLBACK_HANDOVERS if handover isa Integer || h < handover]]
    x = initial_state(p)
    for (k, h) in enumerate(stages)
        last = k == length(stages)
        plan_phases!(pr, [h isa Integer ? ("Cycles", NaN, "cycle") : ("Anderson", h, "cycle"), ("Newton", tol, "iteration")])
        loose = h isa Integer ? _plain_cycles(model, x, prm, forces, dt, h, pr; info_level = info_level, kwarg...) :
            simulate_to_cyclic_steady_state(model, x, prm, forces, dt; tol = h, max_cycles = max_cycles,
                progress = pr, info_level = info_level, kwarg...)
        append!(history, loose.history)
        cycles += loose.cycles
        r = try
            newton_cyclic_steady_state(model, loose.state0, prm, forces, dt; tol = tol, max_iterations = 30,
                krylov = krylov, progress = pr, info_level = info_level, kwarg...)
        catch e
            (_is_solver_failure(e) && !last) || rethrow()
            @debug "Newton's method failed" exception = (e, catch_backtrace())
            nothing
        end
        if !isnothing(r)
            append!(history, r.history)
            cycles += r.cycles
        end
        (last || r.converged) && return (result = r, history = history, cycles = cycles)
        x = isnothing(r) || isempty(loose.history) || r.history[end] >= loose.history[end] ? loose.state0 : r.state0
        report_note!(pr, "  Newton's method " * (isnothing(r) ? "failed" : "did not converge") *
            " from " * _handover_text(h) * "; Anderson acceleration on to $(_sci(stages[k + 1])), then Newton's method again")
    end
end

# The same process with half as many cells (rounded up) in every bed, if any
# bed has more than `coarse_cells`; nothing if not, or if the process cannot
# be built on fewer cells (a bed added without add_bed!, or a parameter set
# per cell). Settings made through the constructor, the design variables, the
# feed or the steps carry over; the coarse grid only gives a start, so a
# setting that does not changes the time taken, not the steady state.
function _coarse_copy(p::Process, coarse_cells)
    (isnothing(coarse_cells) || coarse_cells === false) && return nothing
    fs = p.flowsheet
    bs = beds(fs)
    (isempty(bs) || !all(b -> haskey(fs.specs, b), bs)) && return nothing
    maximum(b -> fs.specs[b].ncells, bs) > coarse_cells || return nothing
    coarse = deepcopy(fs)
    for b in bs
        spec = coarse.specs[b]
        spec.ncells = max(cld(spec.ncells, 2), 2)
        network(coarse).units[b] = fixed_bed_model(coarse.gas, coarse.feed, spec.adsorbent, spec.column; ncells = spec.ncells,
            thermal = spec.thermal, dispersion_velocity = spec.velocity_used, dispersion = spec.dispersion,
            permeability = spec.permeability)
    end
    o = p.options
    q = try
        Process(coarse, p.cycle; name = p.name, timestepping = p.timestepping, key_component = o.key_component,
            ambient_pressure = o.ambient_pressure, pump_efficiency = o.pump_efficiency, heat_capacity_ratio = o.heat_capacity_ratio,
            ramp_rate = o.ramp_rate, initial_state = p.initial_conditions, parameters = p.parameter_settings)
    catch e
        e isa Union{MoccaInputError, DimensionMismatch} || rethrow()
        @debug "No coarser copy of the process" exception = (e, catch_backtrace())
        return nothing
    end
    # The same time steps: the same design variables decide whether a
    # duration stretches the pattern of the steps
    append!(q.variables, p.variables)
    q.reference_durations = p.reference_durations
    return q
end

# The coarser grids _newton_steady_state goes through, coarsest first
function _coarse_grids(p::Process, coarse_cells)
    grids = Process[]
    q = _coarse_copy(p, coarse_cells)
    while !isnothing(q)
        pushfirst!(grids, q)
        q = _coarse_copy(q, coarse_cells)
    end
    return grids
end

function _cells_text(p::Process)
    n = [p.flowsheet.specs[b].ncells for b in beds(p.flowsheet) if haskey(p.flowsheet.specs, b)]
    isempty(n) && return "the grid"
    allequal(n) && return "$(first(n)) cells" * (length(n) > 1 ? " per bed" : "")
    return join(["$b $(p.flowsheet.specs[b].ncells)" for b in beds(p.flowsheet)], ", ") * " cells"
end

# The state `x` of the same process on another grid (see _coarse_copy), as a
# start for `p`: each variable of a unit with holdup interpolated linearly,
# cell centre to cell centre, along its uniform cells (mole fractions stay
# normalised, as their weights add up to one); a unit on the same cells is
# taken as it is
function _interpolate_state(p::Process, x)
    out = initial_state(p)
    for u in holdup_units(p.flowsheet), k in keys(p.model[u].primary_variables)
        xc, xf = x[u][k], out[u][k]
        if size(xc) == size(xf)
            xf .= xc
            continue
        end
        nc, nf = size(xc, ndims(xc)), size(xf, ndims(xf))
        for i in 1:nf
            s = (i - 0.5) / nf * nc + 0.5
            j = clamp(floor(Int, s), 1, nc - 1)
            w = clamp(s - j, 0.0, 1.0)
            if xf isa AbstractMatrix
                xf[:, i] .= (1 - w) .* xc[:, j] .+ w .* xc[:, j + 1]
            else
                xf[i] = (1 - w) * xc[j] + w * xc[j + 1]
            end
        end
    end
    return restart_state(p.model, out)
end

_handover_text(h::Integer) = "$h plain $(_plural("cycle", h))"
_handover_text(h) = "Anderson acceleration to $(_sci(h))"

# `n` cycles from `x`, each from the end of the last, to start Newton's method
# from: the fast modes have died out after a few, and Newton's method does the
# rest
function _plain_cycles(model, x, prm, forces, dt, n, pr; info_level = -1, kwarg...)
    start_phase!(pr, "Cycles"; tol = NaN, estimate = false)
    history = Float64[]
    for k in 1:n
        states, ts = simulate_process(MoccaCase(model, dt, forces; state0 = x, parameters = prm);
            info_level = info_level, output_substates = true, kwarg...)
        (!isempty(ts) && isapprox(sum(ts), sum(dt))) || throw(SimulationFailure("simulation of cycle $k failed"))
        x_end = restart_state(model, states[end])
        push!(history, cycle_change(model, x, x_end))
        report_cycle!(pr, k, history[end]; fraction = k / n)
        x = x_end
    end
    finish_phase!(pr, n, true)
    return (state0 = x, history = history, cycles = n)
end

# ------------------------------------------------------------------------------
# Performance
# ------------------------------------------------------------------------------

_last_cycle(r::ProcessResult) = r.cycles == 1 ? eachindex(r.timesteps) : cycle_window(r.timesteps, cycle_time(r.process), r.cycles)

Operation.adsorbent_volume(p::Process) = sum(b -> adsorbent_volume(p.flowsheet[b], p.parameters[b]), beds(p.flowsheet))
adsorbent_mass(p::Process) = sum(b -> sorbent_mass(p.flowsheet[b], p.parameters[b]), beds(p.flowsheet))

"""
    performance(r::ProcessResult) → NamedTuple

Performance of the process over the last cycle of `r`, from the streams
through its devices (by their roles; see [`add_device!`](@ref)):

- `purity`: mole fraction of the key component in the extract;
- `recovery`: fraction of the key component fed that leaves in the extract;
- `specific_energy` [kWh per tonne of key component in the extract]: vacuum
  pumps on every outlet device whose gas leaves below ambient pressure, and
  compressors on every inlet device whose gas enters a bed above it, as in
  Haghpanah et al. (2013);
- `productivity` [mol of key component per m³ of adsorbent per s];
- `energy` [J per cycle], `cycle_time` [s], the component totals [mol per
  cycle] of the `feed`, `extract` and light `product`, and the
  `adsorbent_volume` [m³] and `adsorbent_mass` [kg].
"""
function performance(r::ProcessResult)
    p = r.process
    fs, o = p.flowsheet, p.options
    w = _last_cycle(r)
    states, ts = r.states, r.timesteps
    total(role) = (ds = devices_with_role(fs, role); N = number_of_components(fs.gas);
        isempty(ds) ? zeros(N) : sum(d -> stream_totals(states, ts, d; window = w), ds))
    n_feed, n_extract, n_product = total(:feed), total(:extract), total(:product)
    E = 0.0
    for role in OUTLET_ROLES, d in devices_with_role(fs, role)
        E += vacuum_pump_energy(states, ts, d; discharge_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio, window = w)
    end
    for role in INLET_ROLES, d in devices_with_role(fs, role)
        E += compressor_energy(states, ts, d; suction_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio, window = w)
    end
    k = o.key_component
    T = sum(ts[w])
    V = adsorbent_volume(p)
    return (purity = purity(n_extract; component = k), recovery = recovery(n_extract, n_feed; component = k),
        specific_energy = specific_energy_kwh_per_tonne(E, n_extract, fs.gas.molar_masses[k]; component = k),
        productivity = productivity(n_extract, T, V; component = k), energy = E, cycle_time = T,
        feed = n_feed, extract = n_extract, product = n_product, adsorbent_volume = V, adsorbent_mass = adsorbent_mass(p))
end

"""
    performance(p::Process; kwarg...) → NamedTuple

Performance at cyclic steady state; keywords go to
[`cyclic_steady_state`](@ref).
"""
performance(p::Process; kwarg...) = performance(cyclic_steady_state(p; kwarg...))

# ------------------------------------------------------------------------------
# Inspection
# ------------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", p::Process)
    fs = p.flowsheet
    println(io, "Process", isempty(p.name) ? "" : " $(p.name)", ": ", join(fs.gas.names, "/"), " feed at ",
        _bar(fs.feed.pressure), ", ", fs.feed.temperature, " K, capturing ", fs.gas.names[p.options.key_component])
    print_flowsheet_diagram(io, fs)
    sts = try
        stages(p)
    catch e
        println(io, "  cycle: ", sprint(showerror, e))
        return
    end
    println(io, "  cycle of $(length(sts)) stages, $(round(sum(s -> s.duration, sts), sigdigits = 5)) s:")
    for st in sts
        open_ = [string(d) for d in unit_names(fs) if haskey(st.settings, d)]
        println(io, "    ", rpad(st.name, 34), lpad(string(round(st.duration, sigdigits = 5)), 8), " s   ", isempty(open_) ? "all closed" : join(open_, ", "))
    end
    if !isempty(p.variables)
        print(io, "  variables: ", join([string(v.name) for v in p.variables], ", "))
    end
end
