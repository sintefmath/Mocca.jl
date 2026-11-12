# Optimisation of a process over its design variables: with gradients, BFGS in
# the box of the variables' bounds with the targets on measures met through an
# augmented Lagrangian; Pareto fronts by sweeping a target; and, where the
# gradients do not help, derivative-free search and exploration of the box.

# The free variables of an optimisation and the map between the unit box of
# their bounds and their values
struct _Design
    names::Vector{Symbol}
    lo::Vector{Float64}
    hi::Vector{Float64}
    x0::Vector{Float64}
    free::Vector{Int}
end

function _design(p::Process, variables, caller)
    names = collect(Symbol, variables)
    vars = [variable(p, n) for n in names]
    for v in vars
        isfinite(v.lower) && isfinite(v.upper) || input_error("variable $(v.name)", "bounds", "$caller needs finite bounds")
    end
    free = findall(v -> v.lower < v.upper, vars)
    isempty(free) && input_error(caller, "variables", "no variable is free to change")
    lo = [v.lower for v in vars]
    hi = [v.upper for v in vars]
    for (k, v) in enumerate(vars)
        k in free || (p[v.name] = v.lower)
    end
    return _Design(names, lo, hi, clamp.(variable_values(p, names), lo, hi), free)
end

# Clamped, as the optimiser can step a round-off outside the unit box
function _to_x(d::_Design, u)
    x = copy(d.x0)
    x[d.free] .= clamp.(d.lo[d.free] .+ u .* (d.hi[d.free] .- d.lo[d.free]), d.lo[d.free], d.hi[d.free])
    return x
end
_to_u(d::_Design, x) = (x[d.free] .- d.lo[d.free]) ./ (d.hi[d.free] .- d.lo[d.free])

# Linear constraints `(Dict(name => coefficient, ...), bound)`, meaning
# Σ coefficient·variable ≤ bound, as A u ≤ b in the unit box
function _linear_constraints(d::_Design, linear_constraints, caller)
    A = zeros(length(linear_constraints), length(d.free))
    b = zeros(length(linear_constraints))
    for (i, (coeffs, bound)) in enumerate(linear_constraints)
        rhs = Float64(bound)
        for (n, a) in pairs(coeffs)
            j = findfirst(==(Symbol(n)), d.names)
            isnothing(j) && input_error(caller, "linear_constraints", "$n is not among the variables")
            rhs -= a * d.lo[j]
            jf = findfirst(==(j), d.free)
            isnothing(jf) || (A[i, jf] = a * (d.hi[j] - d.lo[j]))
        end
        b[i] = rhs
    end
    return (A = A, b = b)
end

# The targets on measures, each a constraint c ≥ 0 in units of the measure's
# scale: c = (m − target)/scale for `at_least`, (target − m)/scale for `at_most`
function _targets(at_least, at_most, sc)
    t = NamedTuple{(:measure, :sign, :target, :scale), Tuple{Symbol, Float64, Float64, Float64}}[]
    for (k, v) in pairs(at_least)
        push!(t, (measure = k, sign = 1.0, target = Float64(v), scale = sc(k)))
    end
    for (k, v) in pairs(at_most)
        push!(t, (measure = k, sign = -1.0, target = Float64(v), scale = sc(k)))
    end
    return t
end
_constraint(t, metrics) = t.sign * (metrics[t.measure] - t.target) / t.scale
_violation(targets, metrics) = maximum((max(0.0, -_constraint(t, metrics)) for t in targets); init = 0.0)

# `f()`, or nothing if it fails because the design cannot be built, such as a
# step that should lower the pressure but would raise it, or because the
# process fails to simulate there
function _try_design(f, values)
    try
        return f()
    catch e
        (e isa MoccaInputError || e isa SimulationFailure || Framework._is_solver_failure(e)) || rethrow()
        @warn "No result at $values" exception = e
        return nothing
    end
end

# Evaluations of the process at designs x, each once, under `key` (the point
# of the unit box the optimiser asked for): its performance, and its gradient
# if asked for. A design whose steady state does not settle, or whose
# simulation fails, is flagged, with no performance.
function _evaluator(p::Process, d::_Design, targets, history, log; gradients = true)
    cache = Dict{Vector{Float64}, Any}()
    function evaluate(x, key = x)
        haskey(cache, key) && return cache[key]
        r = _try_design(NamedTuple{Tuple(d.names)}(Tuple(x))) do
            set_variables!(p, d.names, x)
            if gradients
                performance_gradient(p; variables = d.names, progress = false)
            else
                ss = steady_state(p; progress = false)
                (metrics = performance(ss), gradients = nothing, steady_state = ss)
            end
        end
        settled = !isnothing(r) && r.steady_state.converged === true
        vals = NamedTuple{Tuple(d.names)}(Tuple(x))
        isnothing(r) || settled || @warn "The process did not settle at $vals; the design is left out"
        e = (values = vals, key = copy(key), metrics = isnothing(r) ? nothing : r.metrics, gradients = isnothing(r) ? nothing : r.gradients,
            settled = settled, violation = settled ? _violation(targets, r.metrics) : Inf)
        push!(history, e)
        isnothing(log) || isnothing(r) || log(vals, r)
        cache[copy(key)] = e
        return e
    end
    return evaluate
end

"""
    optimise!(p; variables = all, maximise = nothing, minimise = nothing,
        at_least = (;), at_most = (;), linear_constraints = [], scale = (;),
        method = :bfgs, constraints = :lagrangian, constraint_tol = 1e-4, weight = 1e3,
        max_it = 25, outer_iterations = 5, log = nothing, print_level = 0) → NamedTuple

Optimise the process at the state it settles to: maximise or minimise one
[`performance`](@ref) measure, such as `maximise = :recovery`, subject to
others held above or below targets, such as `at_least = (purity = 0.9,)`, in
the box of the variables' bounds. `linear_constraints` are
`(Dict(name => coefficient, ...), bound)` pairs, meaning
`Σ coefficient·variable ≤ bound`, such as `(Dict(:P_L => 1.0, :P_I => -1.0),
-0.01bar)` for P_I ≥ P_L + 0.01 bar. Measures are compared in units of
`scale` (per measure, 1 by default), so the objective should be of order one.

With `method = :bfgs`, the gradients come from [`performance_gradient`](@ref)
and BFGS in the box (`Jutul.unit_box_bfgs`) takes up to `max_it` iterations
at a time. The targets are met with an augmented Lagrangian (`constraints =
:lagrangian`): each BFGS run maximises the objective less a quadratic
penalty of `weight` on the shortfall from the targets, shifted by a
multiplier per target; after each run the multipliers grow by the shortfall
that is left, and the penalty tenfold where it did not shrink enough. Runs
follow one another, each from the best design of the last, until every
target is met to `constraint_tol` (in units of its scale) and a run
converges before `max_it` iterations, or after `outer_iterations` runs. The
first run is the penalty method, `constraints = :penalty`, which stops there
and may miss the targets slightly.

With `method = :pattern`, the search needs no gradients: a compass search in
the box takes steps along each variable in turn, halving them when none
improves, with a penalty of `weight` on the shortfall from each target, up to
`max_it` evaluations. It suits objectives with kinks or jumps, such as near a
change of regime, where the gradients mislead.

Variables whose bounds are equal stay fixed. Each evaluation finds the steady
state from the last one. A design at which the process does not settle,
fails to simulate or cannot be built (such as an evacuation that would raise
the pressure) is left out, and the optimiser steps back from it.
`log(values, result)` is called after each evaluation. Returns the best
`values` (a NamedTuple), its `metrics`, whether it is `feasible`, the
`multipliers` of the targets, the number of `evaluations` and the `history`
of every evaluation (`values`, `metrics`, `settled` and `violation`, the
largest shortfall from a target in units of its scale), and leaves the
process set to the best values: the best objective among the designs that
meet the targets, or the design nearest to meeting them.
"""
function optimise!(p::Process; variables = variable_names(p), maximise = nothing, minimise = nothing,
        at_least = NamedTuple(), at_most = NamedTuple(), linear_constraints = [], scale = NamedTuple(),
        method = :bfgs, constraints = :lagrangian, constraint_tol = 1e-4, weight = 1e3, max_it = 25,
        outer_iterations = 5, log = nothing, print_level = 0)
    count(!isnothing, (maximise, minimise)) == 1 || input_error("optimise!", "", "give one of maximise or minimise")
    method in (:bfgs, :pattern) || input_error("optimise!", "method", "use :bfgs or :pattern, got $(repr(method))")
    constraints in (:lagrangian, :penalty) || input_error("optimise!", "constraints", "use :lagrangian or :penalty, got $(repr(constraints))")
    objective, sgn = isnothing(maximise) ? (minimise, -1.0) : (maximise, 1.0)
    d = _design(p, variables, "optimise!")
    reset_step_pattern!(p)
    sc(k) = Float64(get(scale, k, 1.0))
    targets = _targets(at_least, at_most, sc)
    lin = _linear_constraints(d, linear_constraints, "optimise!")
    history = Any[]
    evaluate = _evaluator(p, d, targets, history, log; gradients = method == :bfgs)
    F(e) = sgn * e.metrics[objective] / sc(objective)
    λ = zeros(length(targets))
    if method == :pattern
        merit_pattern(u) = (e = evaluate(_to_x(d, u), u); e.settled ? F(e) - weight * sum((max(0.0, -_constraint(t, e.metrics)) for t in targets); init = 0.0) : NaN)
        feasible(u) = isempty(lin.b) || all(lin.A * u .<= lin.b .+ 1e-12)
        _compass_search(merit_pattern, _to_u(d, d.x0), feasible; max_evaluations = max_it)
    else
        ρ = 2weight
        # The augmented Lagrangian of the targets, c ≥ 0, to maximise:
        # F − Σ (max(0, λ − ρc)² − λ²)/2ρ, and its gradient in the unit box
        function lagrangian(e)
            e.settled || return (NaN, zeros(length(d.free)))
            grad(k) = [e.gradients[k][n] for n in d.names]
            v = F(e)
            dv = sgn .* grad(objective) ./ sc(objective)
            for (j, t) in enumerate(targets)
                c = _constraint(t, e.metrics)
                s = max(0.0, λ[j] - ρ * c)
                v -= (s^2 - λ[j]^2) / 2ρ
                s > 0 && (dv = dv .+ s .* t.sign .* grad(t.measure) ./ t.scale)
            end
            return (v, dv[d.free] .* (d.hi[d.free] .- d.lo[d.free]))
        end
        merit(u) = lagrangian(evaluate(_to_x(d, u), u))
        kw = isempty(linear_constraints) ? Dict{Symbol, Any}() : Dict{Symbol, Any}(:lin_ineq => lin)
        u = _to_u(d, d.x0)
        last_violation = Inf
        for k in 1:(constraints == :penalty ? 1 : outer_iterations)
            first_of_run = length(history) + 1
            _, _, bfgs = Jutul.unit_box_bfgs(u, merit; maximize = true, max_it = max_it, grad_tol = 1e-4, obj_change_tol = 1e-5,
                max_initial_update = 0.05, stepIncreaseTol = 3.0, print = print_level, kw...)
            # Whether the run converged, rather than ran out of iterations
            converged = length(bfgs.val) - 1 < max_it
            # The best design of this run for its Lagrangian, from its start
            # (evaluated in an earlier run) and the designs it tried
            run = [e for e in [evaluate(_to_x(d, u), u); history[first_of_run:end]] if e.settled]
            isempty(run) && break
            best_run = argmax(e -> lagrangian(e)[1], run)
            moved = best_run.key != u
            u = best_run.key
            viol = best_run.violation
            constraints == :penalty && break
            # Done when the targets are met and the run converged, or when a
            # run that met them could not improve
            viol <= constraint_tol && (converged || !moved) && break
            for (j, t) in enumerate(targets)
                λ[j] = max(0.0, λ[j] - ρ * _constraint(t, best_run.metrics))
            end
            viol > constraint_tol && viol > 0.25 * last_violation && (ρ *= 10)
            last_violation = viol
        end
    end
    settled = [e for e in history if e.settled]
    isempty(settled) && throw(SimulationFailure("the process settled at none of the designs tried"))
    ok = [e for e in settled if e.violation <= constraint_tol]
    best = isempty(ok) ? argmin(e -> e.violation, settled) : argmax(F, ok)
    if constraints == :penalty && method == :bfgs
        # The penalty method's best: the best penalised objective
        penalised(e) = F(e) - weight * sum((max(0.0, -_constraint(t, e.metrics))^2 for t in targets); init = 0.0)
        best = argmax(penalised, settled)
    end
    set_variables!(p, d.names, collect(values(best.values)))
    return (values = best.values, metrics = best.metrics, feasible = best.violation <= constraint_tol,
        multipliers = NamedTuple{Tuple(t.measure for t in targets)}(Tuple(λ)), evaluations = length(history),
        history = [(values = e.values, metrics = e.metrics, settled = e.settled, violation = e.violation) for e in history])
end

# Compass search in the unit box for the maximum of f: steps of `step` along
# each coordinate in turn, the first that improves taken, the step halved when
# none does, down to `min_step` or `max_evaluations`. A trial point outside
# `feasible`, or where f is not finite, does not improve.
function _compass_search(f, u0, feasible; step = 0.25, min_step = 1e-3, max_evaluations = 50)
    u = copy(u0)
    v = f(u)
    evaluations = 1
    while step >= min_step && evaluations < max_evaluations
        improved = false
        for i in eachindex(u), s in (1.0, -1.0)
            t = copy(u)
            t[i] = clamp(t[i] + s * step, 0.0, 1.0)
            (t == u || !feasible(t)) && continue
            vt = f(t)
            evaluations += 1
            if isfinite(vt) && (!isfinite(v) || vt > v)
                u, v, improved = t, vt, true
                break
            end
            evaluations >= max_evaluations && break
        end
        improved || (step /= 2)
    end
    return (u = u, value = v, evaluations = evaluations)
end

"""
    pareto_front!(p; maximise = nothing, minimise = nothing, sweep, sense = :at_least,
        at_least = (;), at_most = (;), kwarg...) → Vector{NamedTuple}

A Pareto front of two measures by ε-constraints: the optimum of the
objective (`maximise` or `minimise`) for each target in `sweep`, a
`measure => targets` pair such as `:productivity => [0.3, 0.4, 0.5]`, held
`at_least` or, with `sense = :at_most`, at most. The other targets and
keywords go to [`optimise!`](@ref). Each optimisation starts from the optimum
of the one before, so sweep the targets in the order that tightens them.
Returns, per target, the `target`, the optimum's `values` and `metrics`,
whether it is `feasible`, and the number of `evaluations`, and leaves the
process at the last optimum.
"""
function pareto_front!(p::Process; maximise = nothing, minimise = nothing, sweep, sense = :at_least,
        at_least = NamedTuple(), at_most = NamedTuple(), kwarg...)
    sense in (:at_least, :at_most) || input_error("pareto_front!", "sense", "use :at_least or :at_most, got $(repr(sense))")
    measure, sweep_targets = sweep
    return map(sweep_targets) do t
        c = NamedTuple{(Symbol(measure),)}((Float64(t),))
        r = sense == :at_least ?
            optimise!(p; maximise = maximise, minimise = minimise, at_least = merge(at_least, c), at_most = at_most, kwarg...) :
            optimise!(p; maximise = maximise, minimise = minimise, at_least = at_least, at_most = merge(at_most, c), kwarg...)
        (target = Float64(t), values = r.values, metrics = r.metrics, feasible = r.feasible, evaluations = r.evaluations)
    end
end

"""
    explore(p; variables = all, samples = 10, warm_start = false, log = nothing) → Vector{NamedTuple}

The performance of the process at `samples` designs spread over the box of
the variables' bounds (the points of a Halton sequence), to find the regimes
of the process and a start for [`optimise!`](@ref) before following
gradients. Each design is settled from the initial state, unless
`warm_start`, so that one regime does not carry over to the next design.
Returns, per design, its `values`, `metrics` (`nothing` if it could not be
built or failed to simulate) and whether it `settled`, and leaves the process as it was, with
the steady state it starts the next solve from.
`log(values, metrics)` is called after each.
"""
function explore(p::Process; variables = variable_names(p), samples = 10, warm_start = false, log = nothing)
    names = collect(Symbol, variables)
    start = variable_values(p, names)
    warm = (p.warm.state, p.warm.jacobian)
    d = _design(p, names, "explore")
    points = map(1:samples) do k
        u = [_halton(k, _PRIMES[i]) for i in eachindex(d.free)]
        x = _to_x(d, u)
        vals = NamedTuple{Tuple(names)}(Tuple(x))
        ss = _try_design(vals) do
            set_variables!(p, names, x)
            steady_state(p; warm_start = warm_start, progress = false)
        end
        metrics = isnothing(ss) ? nothing : performance(ss)
        isnothing(log) || log(vals, metrics)
        (values = vals, metrics = metrics, settled = !isnothing(ss) && ss.converged === true)
    end
    set_variables!(p, names, start)
    p.warm.state, p.warm.jacobian = warm
    return points
end

const _PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)

# Point k of the Halton sequence in base b: the digits of k mirrored about the
# radix point
function _halton(k, b)
    x, f = 0.0, 1.0
    while k > 0
        f /= b
        x += f * (k % b)
        k ÷= b
    end
    return x
end
