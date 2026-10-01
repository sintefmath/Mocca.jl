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
    CompressorObjective(device; suction_pressure = 1e5, efficiency = 0.72, γ = 1.4, scale = 1.0)

Objective for adjoint gradients: `scale` times the work [J] to compress the
gas entering through `device` from `suction_pressure` up to the pressure at
its outlet, such as a feed compressor taking flue gas at 1 bar up to the bed
inlet pressure. Only steps with forward flow and an outlet above the suction
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
CompressorObjective(device::Symbol; suction_pressure = 1e5, efficiency = 0.72, γ = 1.4, scale = 1.0) =
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
    (!isempty(dt) && isapprox(sum(dt), sum(timesteps))) || error("The simulation failed at t = $(sum(dt)) s")
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
# `time_scale` force), and the rate of each pressure ramp, which runs on time
# since the start of its stage. Later stages run on their own time, so they
# are unchanged. The objective's own dependence on the time step is added
# separately.

_scale_time(x, s) = x
_scale_time(r::ExponentialRamp, s) = ExponentialRamp(r.start, r.stop, r.λ * s; t0 = r.t0)
_scale_time(f::PortSetpoint{K}, s) where K =
    PortSetpoint{K}(PortCondition(_scale_time(f.condition.pressure, s), f.condition.temperature, f.condition.composition))
function _scale_time(f::NamedTuple, s)
    g = map(v -> _scale_time(v, s), f)
    return haskey(g, :time_scale) ? merge(g, (time_scale = s,)) : g
end
_scale_time(f::AbstractDict, s) = Dict{Symbol, Any}(k => _scale_time(v, s) for (k, v) in f)

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
    function F(X, step_info)
        s = X[n_forces .+ (1:length(sets))]
        forces_X = [_scale_time(devectorize(X, i), s[i]) for i in eachindex(sets)]
        return JutulCase(model, dt, forces_X[set_of_step]; state0 = state0, parameters = parameters)
    end
    # One setup serves every objective: the residual's dependence on X, and
    # so its sparsity, does not depend on the objective. Sparsity from every
    # step, so that every stage's forces are seen.
    objectives = G isa NamedTuple ? G : (G = G,)
    packed = Jutul.AdjointPackedResult(states, dt, forces)
    storage = Jutul.AdjointsDI.setup_adjoint_storage_generic(X, F, packed, first(objectives);
        state0 = state0, info_level = -1, single_step_sparsity = false, sparsity_step_type = :all)
    N = length(dt)
    t = cumsum(dt)
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

Extra keywords go to [`setup_process_simulator`](@ref).
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

"""
    newton_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-10, max_iterations = 10, info_level = -1, kwarg...)

Converge a cyclic steady state with Newton's method, from a `state0` near it.
The cycle, given by its `forces` and `timesteps`, maps its initial state `x` to
its final state `Φ(x)`, and each iteration solves

    (I − ∂Φ/∂x) δ = Φ(x) − x,   x ← x + δ.

The Jacobian `∂Φ/∂x` is found row by row with adjoints over the cycle, one for
each primary variable of the units with holdup. That is costly for large
models, but convergence is quadratic, so a few iterations reach a steady state
far tighter than [`simulate_to_cyclic_steady_state`](@ref) can, which needs
many cycles when the slowest mode of the cycle decays slowly. Start it from the
`state0` of a loosely converged `simulate_to_cyclic_steady_state` (`tol` of
`1e-3` to `1e-4`), and start each stage with short steps (`first_dt` in
[`setup_schedule`](@ref)) so that the cycle map is smooth; see
[`cyclic_steady_state_gradient`](@ref).

With `jacobian`, the `jacobian` returned by an earlier call (for example for
slightly different settings, as in an optimisation), the iterations start
from it instead of building one, and improve it with a Broyden update from
each iteration's step and change in residual (a quasi-Newton method), each
costing one cycle instead of one adjoint per variable. The Jacobian is
rebuilt whenever an iteration reduces the cycle change by less than a factor
`1/refresh_ratio`.

Stops when the [`cycle_change`](@ref) is below `tol`. Returns the same fields
as `simulate_to_cyclic_steady_state`, with `cycles` counting the iterations,
plus `jacobian`, the last Jacobian used (as the inverse of `I − ∂Φ/∂x`), and
`jacobian_builds`, how many were built. Extra keywords go to
[`setup_process_simulator`](@ref).
"""
function newton_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-10, max_iterations = 10, jacobian = nothing, refresh_ratio = 0.5, info_level = -1, kwarg...)
    x = restart_state(model, state0)
    history = Float64[]
    states = dt = nothing
    J = isnothing(jacobian) ? nothing : (M = copy(jacobian.M), dofs = jacobian.dofs, map = jacobian.map)
    builds = 0
    previous = nothing
    for it in 0:max_iterations
        states, dt, step_forces = _simulate_steps(model, x, parameters, forces, timesteps; kwarg...)
        x_end = restart_state(model, states[end])
        push!(history, cycle_change(model, x, x_end))
        info_level >= 0 && println("Newton iteration $it: cycle change $(history[end])")
        if history[end] < tol
            return (states = states, timesteps = dt, state0 = x, history = history, cycles = it, converged = true,
                jacobian = J, jacobian_builds = builds)
        end
        it == max_iterations && break
        stalled = it > 0 && history[end] > refresh_ratio * history[end - 1]
        if isnothing(J) || stalled
            J = _cycle_jacobian(model, x, parameters, states, dt, step_forces)
            builds += 1
            previous = nothing
        end
        z0 = Jutul.vectorize_variables(model, x, J.map)
        z1 = Jutul.vectorize_variables(model, x_end, J.map)
        r = z1[J.dofs] .- z0[J.dofs]
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
        # The devices start from the end of the cycle
        z = copy(z1)
        z[J.dofs] .= z0[J.dofs] .+ J.M * r
        x = restart_state(model, _clamp_to_bounds!(model, Jutul.devectorize_variables!(deepcopy(x), model, z, J.map)))
    end
    return (states = states, timesteps = dt, state0 = x, history = history, cycles = max_iterations, converged = false,
        jacobian = J, jacobian_builds = builds)
end

# Jacobian ∂Φ/∂x of the cycle with respect to the initial state of the units
# with holdup, row by row from adjoints, returned as M = (I − ∂Φ/∂x)⁻¹
function _cycle_jacobian(model, x, parameters, states, dt, forces)
    sensitivity, map, = _state_adjoint(model, x, parameters, states, dt, forces)
    dofs = _holdup_dofs(model, map)
    n = length(dofs)
    A = zeros(n, n)
    w = zeros(length(Jutul.vectorize_variables(model, x, map)))
    for (r, i) in enumerate(dofs)
        w .= 0
        w[i] = 1
        A[r, :] .= sensitivity(_CycleObjective(nothing, copy(w), map))[dofs]
    end
    return (M = inv(Matrix{Float64}(I, n, n) .- A), dofs = dofs, map = map)
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

Extra keywords go to [`setup_process_simulator`](@ref).
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
