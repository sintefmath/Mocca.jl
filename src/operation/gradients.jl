# Adjoint gradients of process objectives: with respect to the device settings
# of each stage, and at cyclic steady state.
#
# Objectives are Jutul sum objectives, called for every step as
# G(model, state, dt, step_info, forces) and summed. StreamObjective and
# VacuumPumpObjective give the step terms of stream_totals and
# vacuum_pump_energy; ratios such as purity, recovery and specific energy
# follow from their gradients by the quotient rule.

"""
    StreamObjective(device; component = 1, scale = 1.0)

Objective for adjoint gradients: `scale` times the moles of `component` that
cross `device` from its inlet to its outlet, summed over the steps as in
[`stream_totals`](@ref).
"""
struct StreamObjective
    device::Symbol
    component::Int
    scale::Float64
end
StreamObjective(device::Symbol; component = 1, scale = 1.0) = StreamObjective(device, component, scale)

function (G::StreamObjective)(model, state, dt, step_info, forces)
    s = state[G.device]
    F = s[:MolarFlow][1]
    # Upstream composition. Both sides are evaluated, so the objective reads the
    # same variables whichever way the gas flows.
    y = ifelse(F >= 0, s[:InletComposition][G.component, 1], s[:OutletComposition][G.component, 1])
    return G.scale * dt * F * y
end

"""
    VacuumPumpObjective(device; discharge_pressure = 101325.0, efficiency = 0.72,
        γ = 1.4, scale = 1.0)

Objective for adjoint gradients: `scale` times the vacuum pump work [J] for the
gas leaving through `device`, as in [`vacuum_pump_energy`](@ref).
"""
struct VacuumPumpObjective
    device::Symbol
    discharge_pressure::Float64
    efficiency::Float64
    γ::Float64
    scale::Float64
end
VacuumPumpObjective(device::Symbol; discharge_pressure = 101325.0, efficiency = 0.72, γ = 1.4, scale = 1.0) =
    VacuumPumpObjective(device, discharge_pressure, efficiency, γ, scale)

function (G::VacuumPumpObjective)(model, state, dt, step_info, forces)
    s = state[G.device]
    F = s[:MolarFlow][1]
    P = s[:OutletPressure][1]
    T = s[:InletTemperature][1]
    k = (G.γ - 1) / G.γ
    W = F * GAS_CONSTANT * T / k * ((G.discharge_pressure / P)^k - 1) / G.efficiency
    return G.scale * dt * ifelse((F > 0) & (P < G.discharge_pressure), W, zero(W))
end

"""
    CompressorObjective(device; suction_pressure = 101325.0, efficiency = 0.72, γ = 1.4, scale = 1.0)

Objective for adjoint gradients: `scale` times the work [J] to compress the
gas entering through `device` from `suction_pressure` up to the pressure at
its outlet, such as a feed compressor taking flue gas at atmospheric pressure
up to the bed inlet pressure, as in [`compressor_energy`](@ref). Only steps with forward flow and an outlet above the suction
pressure count. The gas is compressed adiabatically with isentropic
`efficiency`, from the device's inlet temperature.
"""
struct CompressorObjective
    device::Symbol
    suction_pressure::Float64
    efficiency::Float64
    γ::Float64
    scale::Float64
end
CompressorObjective(device::Symbol; suction_pressure = 101325.0, efficiency = 0.72, γ = 1.4, scale = 1.0) =
    CompressorObjective(device, suction_pressure, efficiency, γ, scale)

function (G::CompressorObjective)(model, state, dt, step_info, forces)
    s = state[G.device]
    F = s[:MolarFlow][1]
    P = s[:OutletPressure][1]
    T = s[:InletTemperature][1]
    k = (G.γ - 1) / G.γ
    W = F * GAS_CONSTANT * T / k * ((P / G.suction_pressure)^k - 1) / G.efficiency
    return G.scale * dt * ifelse((F > 0) & (P > G.suction_pressure), W, zero(W))
end

# ------------------------------------------------------------------------------
# Forward simulation for the adjoints
# ------------------------------------------------------------------------------

# Simulate with every step the solver took stored, and return those steps as
# steps of their own. Jutul's adjoints for forces handle report steps split
# into sub-steps incorrectly (they look up the forces of a sub-step by its
# report step), and its DI adjoints use the initial state for every sub-step
# of the first report step, so the adjoints here only ever see unsplit steps.
function _simulate_steps(model, state0, parameters, forces, timesteps; nonlinear_tolerance = 1e-8, info_level = -1, kwarg...)
    sim, cfg = setup_process_simulator(model, state0, parameters;
        info_level = info_level, nonlinear_tolerance = nonlinear_tolerance, output_substates = true, kwarg...)
    result = Jutul.simulate!(sim, timesteps; forces = forces, config = cfg)
    states, dt, step_ix = Jutul.expand_to_ministeps(result)
    (!isempty(dt) && isapprox(sum(dt), sum(timesteps))) || throw(SimulationFailure("the simulation stopped at t = $(sum(dt)) s of $(sum(timesteps)) s"))
    step_forces = forces isa AbstractVector ? forces[step_ix] : fill(forces, length(dt))
    return (states, dt, step_forces)
end

function _evaluate_objective(G, model, states, dt, forces)
    N = length(dt)
    t = cumsum(dt)
    return sum(G(model, states[i], dt[i], Jutul.optimization_step_info(i, t[i], dt[i]; Nstep = N), forces[i]) for i in 1:N)
end

# ------------------------------------------------------------------------------
# Gradients with respect to forces
# ------------------------------------------------------------------------------

# The distinct force sets, in order of first use, and the set of each step.
# Each stage of a schedule shares one force set between its steps.
function _force_sets(forces)
    sets = Any[]
    set_of_step = Int[]
    for f in forces
        k = findfirst(s -> s === f, sets)
        if isnothing(k)
            push!(sets, f)
            k = length(sets)
        end
        push!(set_of_step, k)
    end
    return (sets, set_of_step)
end

# Stage durations. A stage of duration τ is stretched to sτ by stretching each
# of its steps by s. The adjoints differentiate with respect to s at s = 1,
# with the steps kept as simulated: s multiplies the time step in the time
# derivatives of each unit (the TimeScale parameter, set through the units'
# `time_scale` force), and the rate of each pressure ramp over the stage, from
# the stage's start time `t_start`. A ramp may have started before the stage,
# when a step is split into several stages by parallel_stages, so it is first
# restarted at `t_start` from its value there, which leaves it unchanged as
# simulated, and only the part within the stage is stretched. Later stages run
# on their own time, so they are unchanged. The objective's own dependence on
# the time step is added separately.

_scale_time(x, s, t_start) = x
function _scale_time(r::ExponentialRamp, s, t_start)
    return ExponentialRamp(evaluate_profile(r, t_start), r.stop, r.λ * s; t0 = t_start)
end
_scale_time(f::PortSetpoint{K}, s, t_start) where K =
    PortSetpoint{K}(PortCondition(_scale_time(f.condition.pressure, s, t_start), f.condition.temperature, f.condition.composition))
function _scale_time(f::NamedTuple, s, t_start)
    g = map(v -> _scale_time(v, s, t_start), f)
    return haskey(g, :time_scale) ? merge(g, (time_scale = s,)) : g
end
_scale_time(f::AbstractDict, s, t_start) = Dict{Symbol, Any}(k => _scale_time(v, s, t_start) for (k, v) in f)

function _force_adjoint(model, state0, parameters, states, dt, forces, G)
    sets, set_of_step = _force_sets(forces)
    configs = Any[]
    offsets = [0]
    X = Float64[]
    for f in sets
        x, cfg = Jutul.vectorize_forces(f, model)
        append!(X, x)
        push!(configs, cfg)
        push!(offsets, length(X))
    end
    # Then the relative duration of each stage, 1 as simulated
    n_forces = length(X)
    append!(X, ones(length(sets)))
    devectorize(X, i) = Jutul.devectorize_forces(sets[i], model, X[(offsets[i] + 1):offsets[i + 1]], configs[i])
    t = cumsum(dt)
    t_start = [t[k] - dt[k] for k in (findfirst(==(i), set_of_step) for i in eachindex(sets))]
    function F(X, step_info)
        s = X[n_forces .+ (1:length(sets))]
        forces_X = [_scale_time(devectorize(X, i), s[i], t_start[i]) for i in eachindex(sets)]
        return MoccaCase(model, dt, forces_X[set_of_step]; state0 = state0, parameters = parameters)
    end
    # One setup serves every objective: the residual's dependence on X, and
    # so its sparsity, does not depend on the objective. Sparsity from every
    # step, so that every stage's forces are seen.
    objectives = G isa NamedTuple ? G : (G = G,)
    packed = Jutul.AdjointPackedResult(states, dt, forces)
    storage = Jutul.AdjointsDI.setup_adjoint_storage_generic(X, F, packed, first(objectives);
        state0 = state0, info_level = -1, single_step_sparsity = false, sparsity_step_type = :all)
    N = length(dt)
    durations = [sum(dt[set_of_step .== i]) for i in eachindex(sets)]
    out = map(objectives) do Gi
        storage[:adjoint_objective_helper].G = Jutul.adjoint_wrap_objective(Gi, model)
        dX = zeros(length(X))
        Jutul.AdjointsDI.solve_adjoint_generic!(dX, X, F, storage, packed, Gi; state0 = state0, info_level = -1)
        ds = dX[n_forces .+ (1:length(sets))]
        # The objective's own dependence on the time step: each step of stage
        # i lasts s·dt
        for k in 1:N
            info = Jutul.optimization_step_info(k, t[k], dt[k]; Nstep = N)
            ∂G∂dt = Jutul.ForwardDiff.derivative(h -> Gi(model, states[k], h, info, forces[k]), dt[k])
            ds[set_of_step[k]] += dt[k] * ∂G∂dt
        end
        (forces = [devectorize(dX, i) for i in eachindex(sets)], durations = ds ./ durations)
    end
    return G isa NamedTuple ? out : only(out)
end

"""
    force_gradients(model, state0, parameters, forces, timesteps, G;
        nonlinear_tolerance = 1e-8, kwarg...) → NamedTuple

Gradient of the objective `G` with respect to the device settings and the
duration of each stage, from one simulation from `state0`. `G` is a Jutul sum
objective, such as a [`StreamObjective`](@ref).

Returns the `objective` value and, with one entry per stage (per distinct
force set, in order of use):
- `forces`: each entry has the layout of the stage's forces, with every
  number replaced by the derivative of the objective with respect to it:
  `g.forces[2][:V_feed].law.rate` is the derivative with respect to the feed
  rate in the second stage. Mole fractions are treated as independent.
- `durations`: the derivative of the objective with respect to the duration of
  the stage [per s]. The stage is stretched with its steps in proportion;
  pressure ramps run on the time since the start of their stage, so the later
  stages are the same in their own time. This is the derivative of the
  simulation as discretised, the change seen with `reference_durations` in
  [`setup_schedule`](@ref); without it, a new duration keeps the stage's
  first steps, which differs by the time-stepping error. An objective that
  depends on the time other than through the time step is differentiated as
  if it did not.

Extra keywords go to `setup_process_simulator`.
"""
function force_gradients(model, state0, parameters, forces, timesteps, G; kwarg...)
    states, dt, step_forces = _simulate_steps(model, state0, parameters, forces, timesteps; kwarg...)
    g = _force_adjoint(model, state0, parameters, states, dt, step_forces, G)
    return (objective = _evaluate_objective(G, model, states, dt, step_forces), forces = g.forces, durations = g.durations)
end

# ------------------------------------------------------------------------------
# Gradients at cyclic steady state
# ------------------------------------------------------------------------------

# Objective G plus w · x_end, where x_end is the vector of primary variables at
# the end of the cycle, in Jutul's vectorised order (the order of the initial
# state sensitivities). Either part may be absent.
struct _CycleObjective{F, W, M}
    G::F
    w::W
    map::M
end

function (o::_CycleObjective)(model, state, dt, step_info, forces)
    v = isnothing(o.G) ? 0.0 : o.G(model, state, dt, step_info, forces)
    if !isnothing(o.w) && step_info[:step] == step_info[:Nstep]
        v += _dot_primary(o.w, model, state, o.map)
    end
    return v
end

function _dot_primary(w, model::Jutul.MultiModel, state, map)
    return sum(k -> _dot_primary(w, model[k], state[k], map[k]), Jutul.submodels_symbols(model))
end

function _dot_primary(w, model::Jutul.SimulationModel, state, map)
    v = 0.0
    for (name, info) in map
        x = state[name]
        m = info.n_row
        n = info.n_full ÷ m
        # Entities of each row in turn; a fraction variable has one row fewer
        # than it has components
        for j in 1:m, i in 1:n
            v += w[info.offset_x + (j - 1) * n + i] * (x isa AbstractVector ? x[i] : x[j, i])
        end
    end
    return v
end

# Parameter targets: none, so that an adjoint only gives the initial state
# sensitivity
_no_targets(model::Jutul.MultiModel) = Dict(k => Symbol[] for k in Jutul.submodels_symbols(model))
_no_targets(model::Jutul.SimulationModel) = Symbol[]

# Adjoint for initial state sensitivities over the cycle simulated in `states`.
# Objective sparsity is found once and kept, but the objective here changes
# with w, so the objective gradients are taken dense.
function _state_adjoint(model, x0, parameters, states, dt, forces)
    storage = Jutul.setup_adjoint_storage(model; state0 = x0, parameters = parameters,
        include_state0 = true, use_sparsity = false, targets = _no_targets(model))
    function sensitivity(Gc)
        Jutul.solve_adjoint_sensitivities(model, states, dt, Gc;
            storage = storage, state0 = x0, forces = forces, raw_output = true, info_level = -1)
        return copy(storage.dstate0)
    end
    return (sensitivity, storage.state0_map, storage.backward.model)
end

# Positions in the vectorised primary variables of the units with holdup.
# Device variables are copies of their ports and the flow; the cycle does not
# depend on their initial values.
function _holdup_dofs(model::Jutul.MultiModel, map)
    dofs = Int[]
    for (k, m) in _holdup_units(model), info in values(map[k])
        append!(dofs, info.offset_x .+ (1:info.n_x))
    end
    return dofs
end
_holdup_dofs(model::Jutul.SimulationModel, map) = collect(1:sum(info -> info.n_x, values(map)))

# Forward (tangent) sensitivities of the cycle simulated in `states`: the
# products ∂Φ/∂x · V of the cycle map with the columns of V, in Jutul's
# vectorised order of the primary variables (that of the initial state
# sensitivities). The change of the start carries through each step n as
#
#     ∂Rₙ/∂xₙ dxₙ = −∂Rₙ/∂xₙ₋₁ dxₙ₋₁,
#
# with the Jacobians of the steps as simulated: those of a forward simulator,
# and of one in reverse mode for the dependence on the previous state. Each
# step is assembled and factorised once, on the first product, and kept, so
# that a block of columns, and every later product, costs only the solves.
# Jutul's linear systems are in units of each variable's scale. (Adjoints give
# products with the transpose instead, one sweep each, which suits a row of
# ∂Φ/∂x or the transposed system of cyclic_steady_state_gradient, but not
# GMRES on the Newton system.)
function _state_tangent(model, x0, parameters, states, dt, forces)
    fsim = Jutul.Simulator(model; state0 = deepcopy(x0), parameters = deepcopy(parameters), mode = :forward, extra_timing = nothing)
    bsim = Jutul.Simulator(model; state0 = deepcopy(x0), parameters = deepcopy(parameters), mode = :reverse, extra_timing = nothing)
    vmap, = Jutul.variable_mapper(model, :primary)
    unit = ones(length(Jutul.vectorize_variables(model, x0, vmap)))
    _fill_variable_scale!(unit, model, vmap)
    packed = Jutul.AdjointPackedResult(states, dt, forces)
    steps = nothing
    function assemble()
        return map(eachindex(dt)) do i
            info = packed[i].step_info
            s0 = i == 1 ? x0 : states[i - 1]
            Jutul.adjoint_reassemble!(fsim, states[i], s0, info[:dt], forces[i], info[:time])
            Jutul.adjoint_reassemble!(bsim, states[i], s0, info[:dt], forces[i], info[:time])
            (lu(fsim.storage.LinearizedSystem.jac), copy(bsim.storage.LinearizedSystem.jac))
        end
    end
    function tangent(V)
        isnothing(steps) && (steps = assemble())
        D = V ./ unit
        for (F, B) in steps
            D = -(F \ (B * D))
        end
        return D .* unit
    end
    return (tangent, vmap)
end

_fill_variable_scale!(s, model::Jutul.MultiModel, map) =
    foreach(k -> _fill_variable_scale!(s, model[k], map[k]), Jutul.submodels_symbols(model))
function _fill_variable_scale!(s, model::Jutul.SimulationModel, map)
    for info in values(map)
        isnothing(info.scale) || (s[info.offset_x .+ (1:info.n_x)] .= info.scale)
    end
    return s
end

# Scale of each of `dofs`, as the cycle change measures it: the largest
# magnitude of its variable in `x`, per component for mole fractions and
# loadings, at least the variable's floor
function _dof_scale(model, x, map, dofs)
    s = deepcopy(x)
    for (k, m) in _holdup_units(model), p in keys(m.primary_variables)
        v = _unit_state(s, k)[p]
        sc = _variable_scale(m, p, v)
        v .= v isa AbstractMatrix ? repeat(sc, 1, size(v, 2)) : sc
    end
    return Jutul.vectorize_variables(model, s, map)[dofs]
end

"""
    newton_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-10, max_iterations = 10, krylov = :auto, growth = 100, max_backtracks = 6,
        keep_steps_below = 1e-3, progress = false, info_level = -1, kwarg...)

Converge a cyclic steady state with Newton's method from `state0`, such as the
end of a few plain cycles.
The cycle, given by its `forces` and `timesteps`, maps its initial state `x` to
its final state `Φ(x)`, and each iteration solves

    (I − ∂Φ/∂x) δ = Φ(x) − x,   x ← x + δ.

Products with `∂Φ/∂x` are forward sensitivities over the cycle as simulated:
a change in the start carried through the linearised steps. Each step is
assembled and factorised once, so that many products cost little more than
one. Convergence is quadratic, so a few iterations reach a steady state far
tighter than [`simulate_to_cyclic_steady_state`](@ref) can, which needs many
cycles when the slowest mode of the cycle decays slowly. A rough start is
enough, once a few plain cycles have let the fastest modes die out: the
control of the steps described below does the rest. Start each stage with short steps (`first_dt` in
[`setup_schedule`](@ref)) so that the cycle map is smooth; see
[`cyclic_steady_state_gradient`](@ref).

With `krylov = false` the Jacobian `∂Φ/∂x` is built in full, one column for
each primary variable of the units with holdup, in one sweep over the cycle. With `jacobian`, the `jacobian` returned by an earlier call (for example for
slightly different settings, as in an optimisation), the iterations start
from it instead of building one, and improve it with a Broyden update from
each iteration's step and change in residual (a quasi-Newton method), each
costing one cycle. The Jacobian is
rebuilt whenever an iteration reduces the cycle change by less than a factor
`1/refresh_ratio`.

With `krylov = true` the Jacobian is never built: each iteration solves for
the step with GMRES (Jacobian-free Newton–Krylov), with one forward
sensitivity per GMRES iteration, to a relative residual that tightens with
the cycle change. GMRES needs about as many products as the cycle has modes
that do not die out within a cycle, so this is cheaper than the full
Jacobian when there are many more variables than slow modes. The products
also give an approximation of the Jacobian on the directions GMRES explored,
which preconditions the later iterations and is returned as `jacobian`; a
`jacobian` given preconditions the first. At most `max_products` products
are taken per iteration (by default, the number of variables). With
`krylov = :auto` (the default), GMRES is used when the units with holdup
have more than 800 primary variables (a bed of about 130 cells with two
components, 100 with three), and the full Jacobian otherwise; `krylov = n`
sets that number of variables to `n` instead.

Far from the steady state the full step can overshoot: the slow modes of the
cycle make it long. Each step is tried, and taken if the cycle simulates from
it and the norm of the residual `Φ(x) − x` (in the units of the cycle change)
is less than `growth` times that of the best iterate so far. A Newton step often raises the residual of a
few cells before the next one converges, so only steps that make it much
worse are cut back: a quasi-Newton step first gets a new Jacobian, then the
step is shortened, from a quadratic fit of the residual along it, or to a
quarter if the cycle failed to simulate, up to `max_backtracks` times. If no
step is accepted the next iteration starts from the end of the cycle, as a
plain cycle would. The cycle simulated for an accepted step is that of the
next iteration, so a full step costs no extra cycle.

Where the solver has to cut a time step can depend on where the cycle
starts, and a cut that comes and goes between nearby starts makes the cycle
map jump by about the time-stepping error, which near the steady state is
far more than the change left (on the wet LPP VSA, a few steps at the start
of pressurisation from vacuum). Once the cycle change is below
`keep_steps_below`, every later cycle therefore takes the steps of that
iteration's cycle, cuts included (and adds any further cut it needs), so
that the map is smooth. The steady state found is that of those steps, which
differ from the steps given only where the solver cut them.

Stops when the [`cycle_change`](@ref) is below `tol`. With `progress`
(`true`, an `IO` or a [`CycleProgress`](@ref)), each iteration is printed with
its cycle change and its step if shortened, and each Jacobian built, GMRES
solve or step rejected. Returns the same fields
as `simulate_to_cyclic_steady_state`, with `cycles` counting the iterations,
plus `jacobian`, the last Jacobian used (as the inverse of `I − ∂Φ/∂x`),
`jacobian_builds`, how many were built, `products`, the number of GMRES
products, `trials`, the number of steps rejected, and `fallbacks`, the
iterations that fell back to a plain cycle. Extra keywords go to
`setup_process_simulator`.
"""
function newton_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-10, max_iterations = 10, jacobian = nothing, refresh_ratio = 0.5, krylov = :auto, max_products = nothing,
        max_backtracks = 6, growth = 100.0, keep_steps_below = 1e-3, progress = false, info_level = -1, kwarg...)
    krylov = _use_krylov(model, krylov)
    pr = cycle_progress(progress; info_level = info_level)
    start_phase!(pr, "Newton"; tol = tol, noun = "iteration", estimate = false)
    x = restart_state(model, state0)
    history = Float64[]
    J = isnothing(jacobian) ? nothing : (M = copy(jacobian.M), dofs = jacobian.dofs, map = jacobian.map)
    builds = products = trials = fallbacks = 0
    previous = nothing
    note() = _jacobian_note(builds, products, jacobian; trials = trials, fallbacks = fallbacks)
    # The steps of the cycles. Where the solver has to cut a step depends on
    # where the cycle starts, so near the steady state a cut that comes and
    # goes between nearby starts makes the cycle map jump, by about the error
    # of the time stepping, far more than the change left. Once the cycle
    # change is below `keep_steps_below`, Newton's method keeps the steps of
    # its cycle, cuts included, for every later cycle (adding any further
    # cut), so that the map it works with is smooth. Further away, the cuts
    # differ from those at the steady state, so the steps given are used.
    sched_forces, sched_dt = forces, timesteps
    sim = _simulate_steps(model, x, parameters, sched_forces, sched_dt; kwarg...)
    step = 1.0
    # Residual Φ(x) − x of the best iterate so far, for the step control
    best = nothing
    for it in 0:max_iterations
        states, dt, step_forces = sim
        x_end = restart_state(model, states[end])
        push!(history, cycle_change(model, x, x_end))
        report_cycle!(pr, it, history[end]; note = step < 1 ? (step == 0 ? "plain cycle" : "step $(_sci(step))") : "")
        if history[end] < keep_steps_below && length(dt) > length(sched_dt)
            sched_dt, sched_forces = dt, step_forces
        end
        if history[end] < tol
            finish_phase!(pr, it, true; note = note())
            return (states = states, timesteps = dt, state0 = x, history = history, cycles = it, converged = true,
                jacobian = J, jacobian_builds = builds, products = products, trials = trials, fallbacks = fallbacks)
        end
        it == max_iterations && break
        fresh = false
        function direction()
            if krylov
                tangent, vmap = _state_tangent(model, x, parameters, states, dt, step_forces)
                dofs = isnothing(J) ? _holdup_dofs(model, vmap) : J.dofs
                z0 = Jutul.vectorize_variables(model, x, vmap)
                r = Jutul.vectorize_variables(model, x_end, vmap)[dofs] .- z0[dofs]
                # Relative residual of the linear solve: enough that the next
                # change is about the square of this one, or a tenth of `tol`
                c = history[end]
                rtol = clamp(max(c, 0.1 * tol / c), 1e-12, 1e-2)
                itmax = something(max_products, length(dofs))
                δ, M, n_products = _krylov_step(tangent, length(z0), isnothing(J) ? nothing : J.M, dofs, r,
                    _dof_scale(model, x, vmap, dofs), rtol, itmax, pr)
                products += n_products
                J = (M = M, dofs = dofs, map = vmap)
                fresh = true
                return δ
            end
            stalled = it > 0 && history[end] > refresh_ratio * history[end - 1]
            if isnothing(J) || stalled || fresh === :rebuild
                J = _cycle_jacobian(model, x, parameters, states, dt, step_forces; progress = pr)
                builds += 1
                previous = nothing
                fresh = true
            end
            z0 = Jutul.vectorize_variables(model, x, J.map)
            r = Jutul.vectorize_variables(model, x_end, J.map)[J.dofs] .- z0[J.dofs]
            if !isnothing(previous)
                # Good Broyden update of M ≈ (I − ∂Φ/∂x)⁻¹ from the last step s and
                # the change y in the residual Φ(x) − x, for which M y ≈ −s
                s_k = z0[J.dofs] .- previous.z
                y_k = r .- previous.r
                sM = J.M' * s_k
                denom = dot(sM, y_k)
                abs(denom) > eps() * norm(sM) * norm(y_k) && (J.M .-= (s_k .+ J.M * y_k) * sM' ./ denom)
            end
            previous = (z = z0[J.dofs], r = r)
            return J.M * r
        end
        δ = direction()
        z0 = Jutul.vectorize_variables(model, x, J.map)
        z1 = Jutul.vectorize_variables(model, x_end, J.map)
        # The residual measured in the units of the cycle change; the Newton
        # step reduces its norm at a unit rate
        scale = _dof_scale(model, x, J.map, J.dofs)
        r = z1[J.dofs] .- z0[J.dofs]
        merit0 = norm(r ./ scale)
        # Steps are measured against the best iterate so far, in this one's
        # scale, so that growth cannot compound from one step to the next
        merit_best = isnothing(best) ? Inf : norm(best ./ scale)
        merit0 <= merit_best && (best = r)
        reference = min(merit0, merit_best)
        function trial(t)
            # The devices start from the end of the cycle
            z = copy(z1)
            z[J.dofs] .= z0[J.dofs] .+ t .* δ
            xt = restart_state(model, _clamp_to_bounds!(model, Jutul.devectorize_variables!(deepcopy(x), model, z, J.map)))
            simt = try
                _simulate_steps(model, xt, parameters, sched_forces, sched_dt; kwarg...)
            catch e
                Framework._is_solver_failure(e) || rethrow()
                return (x = xt, sim = nothing, merit = Inf, change = Inf)
            end
            zt = Jutul.vectorize_variables(model, xt, J.map)
            xe = restart_state(model, simt[1][end])
            ze = Jutul.vectorize_variables(model, xe, J.map)
            return (x = xt, sim = simt, merit = norm((ze[J.dofs] .- zt[J.dofs]) ./ scale), change = cycle_change(model, xt, xe))
        end
        # Line search: take the full step if it reduces the residual enough
        # (Armijo), or reaches `tol`. Otherwise a quasi-Newton step first gets
        # a new Jacobian, then the step is cut, from a quadratic fit of the
        # residual along it, or to a quarter if the cycle failed to simulate
        t = 1.0
        accepted = nothing
        for k in 0:max_backtracks
            a = trial(t)
            @debug "Newton step" it t ratio = a.merit / merit0 change = a.change
            if isfinite(a.merit) && (a.change < tol || a.merit <= growth * reference)
                accepted = a
                break
            end
            trials += 1
            count_cycle!(pr)
            why = isfinite(a.merit) ? "residual $(_sci(a.merit / merit0))× the last" : "the cycle failed to simulate"
            if !krylov && !fresh
                report_status!(pr, "step rejected ($why), rebuilding the Jacobian"; persist = true)
                fresh = :rebuild
                δ = direction()
                continue
            end
            k == max_backtracks && break
            t_next = isfinite(a.merit) ? t^2 * merit0^2 / (a.merit^2 - merit0^2 + 2t * merit0^2) : t / 4
            report_status!(pr, "step $(_sci(t)) rejected ($why)"; persist = true)
            t = clamp(t_next, t / 10, t / 2)
        end
        if isnothing(accepted)
            # Every step tried failed: start the next iteration from the end of
            # this cycle, as a plain cycle would
            fallbacks += 1
            report_status!(pr, "no step accepted; continuing from the end of the cycle"; persist = true)
            x = x_end
            sim = _simulate_steps(model, x, parameters, sched_forces, sched_dt; kwarg...)
            step = 0.0
            previous = nothing
        else
            x, sim, step = accepted.x, accepted.sim, t
        end
    end
    finish_phase!(pr, max_iterations, false; note = note())
    return (states = sim[1], timesteps = sim[2], state0 = x, history = history, cycles = max_iterations, converged = false,
        jacobian = J, jacobian_builds = builds, products = products, trials = trials, fallbacks = fallbacks)
end

function _jacobian_note(builds, products, given; trials = 0, fallbacks = 0)
    if products > 0
        note = "$products GMRES $(_plural("product", products))"
    elseif builds == 0
        note = isnothing(given) ? "" : "Jacobian reused"
    else
        note = "$builds $(builds == 1 ? "Jacobian" : "Jacobians") built" * (isnothing(given) ? "" : " (the one given went stale)")
    end
    trials > 0 && (note *= ", $trials $(_plural("step", trials)) rejected")
    fallbacks > 0 && (note *= ", $fallbacks plain $(_plural("cycle", fallbacks))")
    return note
end

# Holdup variables above which `krylov = :auto` solves Newton's steps by
# GMRES. GMRES needs about as many products as the cycle has slow modes,
# which grow far more slowly than the variables, while the full Jacobian
# costs a column per variable. Below this the two cost about the same per
# iteration, and which needs fewer iterations varies from grid to grid; the
# exact steps of the full Jacobian won on the wet LPP VSA at 80 cells (640
# variables) and GMRES on every grid from 960 variables, in cold solves of
# the four-stage and wet LPP VSAs on 30 to 320 cells
const KRYLOV_VARIABLES = 800

# Whether Newton's method on `model` solves its steps by GMRES: `krylov` if
# true or false; with `:auto`, or a number of variables, whether the units
# with holdup have more primary variables than that (KRYLOV_VARIABLES)
function _use_krylov(model, krylov)
    krylov isa Bool && return krylov
    n = krylov === :auto ? KRYLOV_VARIABLES : krylov
    n isa Integer || throw(ArgumentError("krylov must be true, false, :auto or a number of variables, got $(repr(krylov))"))
    return _holdup_count(model) > n
end

_holdup_count(model) = length(_holdup_dofs(model, first(Jutul.variable_mapper(model, :primary))))

# Jacobian ∂Φ/∂x of the cycle with respect to the initial state of the units
# with holdup, a column per variable from one sweep of forward sensitivities,
# returned as M = (I − ∂Φ/∂x)⁻¹
function _cycle_jacobian(model, x, parameters, states, dt, forces; progress = nothing)
    tangent, map = _state_tangent(model, x, parameters, states, dt, forces)
    dofs = _holdup_dofs(model, map)
    n = length(dofs)
    report_status!(progress, "building the Jacobian from $n forward sensitivities over $(length(dt)) steps"; persist = true)
    E = zeros(length(Jutul.vectorize_variables(model, x, map)), n)
    for (c, i) in enumerate(dofs)
        E[i, c] = 1
    end
    A = tangent(E)[dofs, :]
    return (M = inv(Matrix{Float64}(I, n, n) .- A), dofs = dofs, map = map)
end

# One Newton step δ from (I − ∂Φ/∂x) δ = r on `dofs`, by GMRES with forward
# sensitivities, in variables divided by `scale` so that the residual norm
# weighs them as the cycle change does. `M`, an approximation of
# (I − ∂Φ/∂x)⁻¹, preconditions it. Returns δ, M updated with the products, and
# the number of products.
function _krylov_step(tangent, n_full, M, dofs, r, scale, rtol, itmax, progress)
    n = length(dofs)
    w = zeros(n_full)
    U = Vector{Float64}[]
    JU = Vector{Float64}[]
    b = r ./ scale
    function apply!(y, v)
        w[dofs] .= v .* scale
        jv = tangent(w)[dofs] ./ scale
        push!(U, copy(v))
        push!(JU, jv)
        y .= v .- jv
        return y
    end
    function status(solver)
        ρ = solver.stats.residuals
        report_status!(progress, "GMRES: product $(length(U)), residual $(_sci(ρ[end] / ρ[1])) → $(_sci(rtol))")
        return false
    end
    op = Jutul.LinearOperators.LinearOperator(Float64, n, n, false, false, apply!)
    # In scaled variables, M becomes S⁻¹ M S
    Ms = isnothing(M) ? Matrix{Float64}(I, n, n) : M .* (scale .\ scale')
    precond = isnothing(M) ? (ldiv = false,) :
        (N = Jutul.LinearOperators.LinearOperator(Float64, n, n, false, false, (y, v) -> mul!(y, Ms, v)), ldiv = false)
    δs, stats = Jutul.Krylov.gmres(op, b; rtol = rtol, atol = 0.0, itmax = itmax, memory = itmax, history = true, callback = status, precond...)
    res = stats.residuals[end] / stats.residuals[1]
    report_status!(progress, "GMRES: $(length(U)) $(_plural("product", length(U))), residual reduced to $(_sci(res))" *
        (stats.solved ? "" : " (not converged, $(stats.status))"); persist = true)
    isempty(U) || _multisecant_update!(Ms, reduce(hcat, U), reduce(hcat, JU))
    return (δs .* scale, Ms .* (scale ./ scale'), length(U))
end

# Update M ≈ (I − J)⁻¹ in place so that it maps (I − J) U to U for the columns
# of U, with JU = J·U, and is unchanged on what is orthogonal to M'U: the
# multisecant form of the good Broyden update. From M = I it is the inverse of
# I − J restricted to the span of U, J taken as zero elsewhere.
function _multisecant_update!(M, U, JU)
    # An orthonormal basis of the directions, and J on it
    F = svd(U)
    keep = F.S .> 1e-10 * F.S[1]
    T = F.V[:, keep] ./ F.S[keep]'
    Q, JQ = U * T, JU * T
    Z = Q .- JQ
    MZ = M * Z
    QM = Q' * M
    C = QM * Z
    cond(C) < 1e12 || return M
    M .+= (Q .- MZ) * (C \ QM)
    return M
end

"""
    cyclic_steady_state_gradient(model, state0, parameters, forces, timesteps, G;
        targets = nothing, with_forces = false, rtol = 1e-8, maxiter = 100,
        nonlinear_tolerance = 1e-8, kwarg...) → NamedTuple

Gradient at cyclic steady state of the objective `G` over one cycle, given by
its `forces` and `timesteps`. `state0` is the state at the start of a cycle at
steady state, such as `state0` from [`newton_cyclic_steady_state`](@ref) or
[`simulate_to_cyclic_steady_state`](@ref). `G` is a Jutul sum objective, such
as a [`StreamObjective`](@ref), or a `NamedTuple` of them, which share the
simulation of the cycle.

The cycle maps its initial state `x` to its final state `Φ(x, p)`, and the
steady state `x*` solves `x* = Φ(x*, p)`, so it moves with any parameter or
setting `p`. By the implicit function theorem, the gradient of
`J(p) = G(x*(p), p)` is

    dJ/dp = ∂J/∂p + μᵀ ∂Φ/∂p,   where   (I − ∂Φ/∂x)ᵀ μ = ∂J/∂x.

Each product with `(∂Φ/∂x)ᵀ` is an adjoint solve over the cycle, and GMRES
finds `μ` from them; `rtol` and `maxiter` control it. Pass the `jacobian` from
[`newton_cyclic_steady_state`](@ref) to precondition it, which cuts the
iterations from tens to a few. The gradient is only as
accurate as the steady state: an error `e` in the cycle change gives a
gradient error of about `e/(1 − ρ)`, where `ρ` is the decay factor per cycle
of the slowest mode. [`newton_cyclic_steady_state`](@ref) converges it far
enough.

The cycle map must be smooth for the gradient to mean anything. Abrupt stage
starts taken in one long step can make it jump: the solver then lands on a
different solution of that step for a slightly different start. Start each
stage with short steps (`first_dt` in [`setup_schedule`](@ref)).

Returns a NamedTuple (one per objective, if `G` is a `NamedTuple`) with:
- `objective`: the value of `G` over the cycle;
- `parameters`: the gradient with respect to the parameters of each unit, laid
  out like the output of `Jutul.solve_adjoint_sensitivities`, for the
  parameters in `targets` (a `Dict` of unit => parameter names) or all of them;
- `forces` and `durations`: if `with_forces`, the gradient with respect to
  each stage's device settings and duration, as in [`force_gradients`](@ref).
  A longer stage also lengthens the cycle; divide by the cycle time for rates;
- `state0`: the multiplier `μ`, laid out like `parameters`: the change in the
  objective at steady state per unit of each initial-state variable added at
  the start of a cycle;
- `iterations` and `converged` for GMRES.

Extra keywords go to `setup_process_simulator`.
"""
function cyclic_steady_state_gradient(model, state0, parameters, forces, timesteps, G;
        targets = nothing, with_forces = false, jacobian = nothing, rtol = 1e-8, maxiter = 100, info_level = -1, kwarg...)
    x0 = restart_state(model, state0)
    states, dt, step_forces = _simulate_steps(model, x0, parameters, forces, timesteps; kwarg...)
    sensitivity, vmap, state_model = _state_adjoint(model, x0, parameters, states, dt, step_forces)
    storage_kwarg = isnothing(targets) ? NamedTuple() : (targets = targets,)
    storage = Jutul.setup_adjoint_storage(model; state0 = x0, parameters = parameters, use_sparsity = false, storage_kwarg...)
    n = length(Jutul.vectorize_variables(model, x0, vmap))
    apply!(y, w) = (y .= w .- sensitivity(_CycleObjective(nothing, copy(w), vmap)))
    op = Jutul.LinearOperators.LinearOperator(Float64, n, n, false, false, apply!)
    # A Jacobian of a nearby cycle, such as the last one from
    # newton_cyclic_steady_state, makes a good preconditioner
    if isnothing(jacobian)
        precond = (ldiv = false,)
    else
        function apply_precond!(y, v)
            y .= v
            y[jacobian.dofs] .= transpose(jacobian.M) * v[jacobian.dofs]
            return y
        end
        precond = (N = Jutul.LinearOperators.LinearOperator(Float64, n, n, false, false, apply_precond!), ldiv = false)
    end
    objectives = G isa NamedTuple ? G : (G = G,)
    # (I − ∂Φ/∂x)ᵀ μ = ∂J/∂x for each objective
    solves = map(objectives) do Gi
        b = sensitivity(_CycleObjective(Gi, nothing, vmap))
        μ, stats = Jutul.Krylov.gmres(op, b; rtol = rtol, atol = 0.0, itmax = maxiter, memory = maxiter,
            verbose = info_level > 0 ? 1 : 0, precond...)
        stats.solved || @warn "GMRES for the cyclic steady state gradient did not converge in $maxiter iterations: $(stats.status)"
        (μ = μ, stats = stats, Gμ = _CycleObjective(Gi, μ, vmap))
    end
    # Then ∂J/∂p + μᵀ ∂Φ/∂p: one adjoint for the parameters, and one setup of
    # the force adjoints shared by all objectives
    gf = with_forces ? _force_adjoint(model, x0, parameters, states, dt, step_forces, map(s -> s.Gμ, solves)) :
        map(_ -> (forces = nothing, durations = nothing), objectives)
    out = map(objectives, solves, gf) do Gi, s, f
        ∇p = Jutul.solve_adjoint_sensitivities(model, states, dt, s.Gμ;
            storage = storage, state0 = x0, forces = step_forces, raw_output = true, info_level = -1)
        (objective = _evaluate_objective(Gi, model, states, dt, step_forces),
            parameters = Jutul.store_sensitivities(storage.parameter.model, ∇p, storage.parameter_map),
            forces = f.forces, durations = f.durations, state0 = Jutul.store_sensitivities(state_model, s.μ, vmap),
            iterations = s.stats.niter, converged = s.stats.solved)
    end
    return G isa NamedTuple ? out : only(out)
end
