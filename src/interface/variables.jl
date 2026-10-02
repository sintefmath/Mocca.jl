# Named design variables of a Process, their gradients at cyclic steady state,
# and optimisation over them.
#
# A variable sets one or more places: a field of a step of the cycle, of the
# feed, or a parameter of a unit. The gradient of a performance measure with
# respect to the variable follows by the chain rule: the adjoint gives it with
# respect to each stage's device settings and duration, and the change of
# those with the variable comes from compiling the cycle with the variable
# perturbed, which is cheap and, for the built-in steps, linear.

abstract type Place end

"""
    StepPlace(step, field; sequence = "")

The `field` (such as `:to`, `:velocity` or `:duration`) of a step of the
cycle, picked by its position (counting through all steps of the sequence,
including those in `Concurrent` groups) or its name, such as
`"evacuation"`, or `"A evacuation"` for a step of bed `A`. In a cycle of
several sequences, `sequence` is the label of the sequence.
"""
struct StepPlace <: Place
    sequence::String
    step::Union{Int, String}
    field::Symbol
end
StepPlace(step, field; sequence = "") = StepPlace(String(sequence), step isa Integer ? Int(step) : String(step), Symbol(field))

"""
    FeedPlace(field)

A field of the flowsheet's [`Feed`](@ref): `:pressure` or `:temperature`.
"""
struct FeedPlace <: Place
    field::Symbol
end

"""
    ParameterPlace(unit, name)

The Jutul parameter `name` of `unit`, such as `ParameterPlace(:Bed,
:AdsorbentDensity)`, set to one value throughout the unit.
"""
struct ParameterPlace <: Place
    unit::Symbol
    name::Symbol
end

_place(x::Place) = x
_place(x::Tuple{Symbol, Symbol}) = first(x) == :feed ? FeedPlace(last(x)) : StepPlace(String(first(x)), last(x))
_place(x::Tuple{Any, Any}) = StepPlace(x[1], x[2])
_place(x::Tuple{Symbol, Any, Any}) = x[1] == :parameter ? ParameterPlace(Symbol(x[2]), Symbol(x[3])) : StepPlace(x[2], x[3]; sequence = String(x[1]))
_place(x::Tuple{AbstractString, Any, Any}) = StepPlace(x[2], x[3]; sequence = x[1])
_place(x) = input_error("variable", "places", "a place is a StepPlace, FeedPlace or ParameterPlace, or a tuple such as (\"evacuation\", :to); got $x")

"""
A named design variable: the places it sets, its bounds and a description.
"""
struct DesignVariable
    name::Symbol
    places::Vector{Place}
    lower::Float64
    upper::Float64
    description::String
end

"""
    add_variable!(p, name, places...; bounds = (-Inf, Inf), description = "")

Add a design variable `name` to the process, setting every one of `places`
(see [`StepPlace`](@ref), [`FeedPlace`](@ref), [`ParameterPlace`](@ref), or
tuples such as `("evacuation", :to)`, `(:feed, :pressure)` and
`(:parameter, :Bed, :AdsorbentDensity)`). Since pressures carry over from step
to step, one place is usually enough: the evacuation pressure is also where
the next pressurisation starts. `bounds` are needed by [`optimise!`](@ref).

Read and set it with `p[name]` and `p[name] = value`.
"""
function add_variable!(p::Process, name::Symbol, places...; bounds = (-Inf, Inf), description = "")
    any(v -> v.name == name, p.variables) && input_error("variable $name", "", "the process already has a variable of that name")
    isempty(places) && input_error("variable $name", "places", "give at least one place it sets")
    pl = Place[_place(x) for x in places]
    for x in pl
        v = _at(() -> _get(p, x), "variable $name")
        v isa Real || input_error("variable $name", "places", "$(_describe(x)) is $(repr(v)), not a number; set it to a number first")
    end
    lo, hi = Float64(bounds[1]), Float64(bounds[2])
    lo <= hi || input_error("variable $name", "bounds", "lower bound $lo is above the upper bound $hi")
    push!(p.variables, DesignVariable(name, pl, lo, hi, String(description)))
    return p
end

variable_names(p::Process) = [v.name for v in p.variables]

function variable(p::Process, name::Symbol)
    i = findfirst(v -> v.name == name, p.variables)
    isnothing(i) && input_error("process", "variables", "no variable $name; the process has $(join(variable_names(p), ", "))")
    return p.variables[i]
end

Base.getindex(p::Process, name::Symbol) = _get(p, first(variable(p, name).places))
function Base.setindex!(p::Process, value, name::Symbol)
    v = variable(p, name)
    for x in v.places
        _set!(p, x, value)
    end
    return value
end

"Values of the variables `names`."
variable_values(p::Process, names = variable_names(p)) = [p[n] for n in names]

"Set the variables `names` to `values`."
function set_variables!(p::Process, names, values)
    for (n, x) in zip(names, values)
        p[n] = x
    end
    return p
end

_describe(x::StepPlace) = "field $(x.field) of step $(repr(x.step))" * (isempty(x.sequence) ? "" : " of the $(x.sequence) sequence")
_describe(x::FeedPlace) = "the feed $(x.field)"
_describe(x::ParameterPlace) = "parameter $(x.name) of $(x.unit)"

function _find_step(p::Process, x::StepPlace)
    seqs = p.cycle.sequences
    k = findfirst(s -> first(s) == x.sequence, seqs)
    isnothing(k) && input_error("variable", "places", "the cycle has no sequence $(repr(x.sequence)); it has $(join(repr.(first.(seqs)), ", "))")
    steps = [s for (label, s) in cycle_steps(p.cycle) if label == x.sequence && s isa AbstractStep]
    if x.step isa Int
        1 <= x.step <= length(steps) || input_error("variable", "places", "the cycle has $(length(steps)) steps, so there is no step $(x.step)")
        return steps[x.step]
    end
    names_of(s) = (n = step_name(s); [n; ["$b $n" for b in step_beds(s)]])
    matches = [s for s in steps if x.step in names_of(s)]
    length(matches) == 1 && return matches[1]
    isempty(matches) && input_error("variable", "places", "no step named $(repr(x.step)); steps are $(join(repr.(step_name.(steps)), ", "))")
    input_error("variable", "places", "$(length(matches)) steps are named $(repr(x.step)); pick one by position or by bed, such as \"A $(x.step)\"")
end

function _get(p::Process, x::StepPlace)
    s = _find_step(p, x)
    hasfield(typeof(s), x.field) || input_error("variable", "places", "a $(nameof(typeof(s))) step has no field $(x.field); it has $(join(fieldnames(typeof(s)), ", "))")
    return getfield(s, x.field)
end
function _set!(p::Process, x::StepPlace, v)
    s = _find_step(p, x)
    setfield!(s, x.field, convert(fieldtype(typeof(s), x.field), v))
    p.compiled = nothing
end

function _get(p::Process, x::FeedPlace)
    x.field in (:pressure, :temperature) || input_error("variable", "places", "a feed variable sets :pressure or :temperature, got $(x.field)")
    return getfield(p.flowsheet.feed, x.field)
end
_set!(p::Process, x::FeedPlace, v) = (setfield!(p.flowsheet.feed, x.field, Float64(v)); p.compiled = nothing)

function _get(p::Process, x::ParameterPlace)
    haskey(p.parameters, x.unit) || input_error("variable", "places", "$(x.unit) is not a unit of the flowsheet")
    haskey(p.parameters[x.unit], x.name) || input_error("variable", "places", "$(x.unit) has no parameter $(x.name)")
    return first(p.parameters[x.unit][x.name])
end
function _set!(p::Process, x::ParameterPlace, v)
    p.parameters[x.unit][x.name] .= v
    get!(p.parameter_settings, x.unit, Dict{Symbol, Any}())[x.name] = v
end

# ------------------------------------------------------------------------------
# Gradients
# ------------------------------------------------------------------------------

# A sum of objectives, as one adjoint objective
struct _SumObjective{T}
    terms::T
end
(G::_SumObjective)(model, state, dt, step_info, forces) = sum(g(model, state, dt, step_info, forces) for g in G.terms)

function _performance_objectives(p::Process)
    fs, o = p.flowsheet, p.options
    k = o.key_component
    N = number_of_components(fs.gas)
    extract = devices_with_role(fs, :extract)
    feed = devices_with_role(fs, :feed)
    isempty(extract) && input_error("process", "", "performance needs an extract device")
    isempty(feed) && input_error("process", "", "performance needs a feed device")
    energy = Any[]
    for role in OUTLET_ROLES, d in devices_with_role(fs, role)
        push!(energy, VacuumPumpObjective(d; discharge_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio))
    end
    for role in INLET_ROLES, d in devices_with_role(fs, role)
        push!(energy, CompressorObjective(d; suction_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio))
    end
    return (key_extract = _SumObjective(Tuple(StreamObjective(d; component = k) for d in extract)),
        total_extract = _SumObjective(Tuple(StreamObjective(d; component = i) for d in extract for i in 1:N)),
        key_feed = _SumObjective(Tuple(StreamObjective(d; component = k) for d in feed)),
        energy = _SumObjective(Tuple(energy)))
end

# Stage force vectors and durations of the process as compiled now
function _stage_vectors(p::Process)
    sts = stages(p)
    t = 0.0
    vecs = Vector{Vector{Float64}}()
    for st in sts
        f = _stage_forces(p.flowsheet, p.model, st, t)
        push!(vecs, first(Jutul.vectorize_forces(f, p.model)))
        t += st.duration
    end
    return (vecs, [st.duration for st in sts])
end

# Change of every stage's force vector and duration with each variable, by
# central differences on the compiled cycle
function _design_jacobian(p::Process, vars)
    base_v, base_τ = _stage_vectors(p)
    out = map(vars) do v
        all(pl -> pl isa ParameterPlace, v.places) && return nothing
        x0 = p[v.name]
        h = 1e-6 * (x0 == 0 ? 1.0 : abs(x0))
        p[v.name] = x0 + h
        vp, τp = _stage_vectors(p)
        p[v.name] = x0 - h
        vm, τm = _stage_vectors(p)
        p[v.name] = x0
        (length(vp) == length(vm) == length(base_v) && all(length.(vp) .== length.(base_v)) && all(length.(vm) .== length.(base_v))) ||
            input_error("variable $(v.name)", "", "changing it changes the number of stages of the cycle, so it has no gradient here")
        (dF = [(a .- b) ./ 2h for (a, b) in zip(vp, vm)], dτ = (τp .- τm) ./ 2h)
    end
    _stage_vectors(p)  # leave the process compiled at the base values
    return out
end

"""
    performance_gradient(p; variables = all, steady_state = nothing, rtol = 1e-7) → NamedTuple

The [`performance`](@ref) of the process at cyclic steady state and its
gradient with respect to the design `variables` (names; see
[`add_variable!`](@ref)), by adjoints and the implicit function theorem (see
[`cyclic_steady_state_gradient`](@ref)). Purity, recovery, specific energy and
productivity follow from cycle sums by the quotient rule.

Returns `metrics`, the performance; `gradients`, a NamedTuple per measure of
NamedTuples per variable, such as `g.gradients.recovery.P_L`; and the
`steady_state` used. The axial dispersion and the gas properties of a bed stay
as they were when the process was built.
"""
function performance_gradient(p::Process; variables = variable_names(p), steady_state = nothing, rtol = 1e-7, info_level = -1)
    names = collect(Symbol, variables)
    vars = [variable(p, n) for n in names]
    ss = isnothing(steady_state) ? cyclic_steady_state(p; info_level = info_level) : steady_state
    forces, dt = cycle_schedule(p)
    objectives = _performance_objectives(p)
    need_forces = any(v -> any(pl -> !(pl isa ParameterPlace), v.places), vars)
    targets = _no_targets(p.model)
    for v in vars, pl in v.places
        pl isa ParameterPlace && push!(targets[pl.unit], pl.name)
    end
    g = cyclic_steady_state_gradient(p.model, ss.state0, p.parameters, forces, dt, objectives;
        with_forces = need_forces, jacobian = ss.jacobian, targets = targets, rtol = rtol, info_level = info_level)
    jac = need_forces ? _design_jacobian(p, vars) : fill(nothing, length(vars))
    function d_sum(r)
        map(vars, jac) do v, J
            total = 0.0
            if !isnothing(J)
                for (s, dF) in enumerate(J.dF)
                    gv = first(Jutul.vectorize_forces(r.forces[s], p.model))
                    total += sum(gv .* dF)
                end
                total += sum(r.durations .* J.dτ)
            end
            for pl in v.places
                pl isa ParameterPlace && (total += sum(r.parameters[pl.unit][pl.name]))
            end
            total
        end
    end
    sums = map(r -> r.objective, g)
    dsums = map(d_sum, g)
    dT = [isnothing(J) ? 0.0 : sum(J.dτ) for J in jac]
    a, b, f, E = sums.key_extract, sums.total_extract, sums.key_feed, sums.energy
    da, db, df, dE = dsums.key_extract, dsums.total_extract, dsums.key_feed, dsums.energy
    T = cycle_time(p)
    V = adsorbent_volume(p)
    M = p.flowsheet.gas.molar_masses[p.options.key_component]
    to_kwh_per_t = 1 / 3.6e6 / (M / 1000)
    metrics = performance(ProcessResult(p, ss.states, ss.timesteps, 1, ss.state0, ss.converged, ss.history, ss.iterations, ss.jacobian))
    grads = (purity = (da .* b .- a .* db) ./ b^2, recovery = (da .* f .- a .* df) ./ f^2,
        specific_energy = (dE .* a .- E .* da) ./ a^2 .* to_kwh_per_t,
        productivity = (da .* T .- a .* dT) ./ (V * T^2))
    tup(v) = NamedTuple{Tuple(names)}(Tuple(v))
    return (metrics = metrics, gradients = map(tup, grads), steady_state = ss)
end

# ------------------------------------------------------------------------------
# Optimisation
# ------------------------------------------------------------------------------

"""
    optimise!(p; variables = all, maximise = nothing, minimise = nothing,
        at_least = (;), at_most = (;), linear_constraints = [], scale = (;),
        weight = 1e3, max_it = 25, log = nothing, print_level = 0) → NamedTuple

Optimise the process at cyclic steady state with gradients (BFGS in the
box of the variables' bounds, `Jutul.unit_box_bfgs`): maximise or minimise
one [`performance`](@ref) measure, such as `maximise = :recovery`, subject to
others held above or below targets, such as `at_least = (purity = 0.9,)`.
Measure constraints are quadratic penalties with `weight` on the shortfall in
units of `scale` (per measure, 1 by default), so the result may miss them
slightly. `linear_constraints` are `(Dict(name => coefficient, ...), bound)`
pairs, meaning `Σ coefficient·variable ≤ bound`, such as
`(Dict(:P_L => 1.0, :P_I => -1.0), -0.01bar)` for P_I ≥ P_L + 0.01 bar.

Variables whose bounds are equal stay fixed. Each evaluation finds the cyclic
steady state from the last one. `log(values, result)` is called after each.
Returns the best `values` (a NamedTuple), its `metrics`, the number of
`evaluations` and the `history` of every evaluation, and leaves the process
set to the best values.
"""
function optimise!(p::Process; variables = variable_names(p), maximise = nothing, minimise = nothing,
        at_least = NamedTuple(), at_most = NamedTuple(), linear_constraints = [], scale = NamedTuple(),
        weight = 1e3, max_it = 25, log = nothing, print_level = 0)
    count(!isnothing, (maximise, minimise)) == 1 || input_error("optimise!", "", "give one of maximise or minimise")
    objective, sgn = isnothing(maximise) ? (minimise, -1.0) : (maximise, 1.0)
    names = collect(Symbol, variables)
    vars = [variable(p, n) for n in names]
    for v in vars
        isfinite(v.lower) && isfinite(v.upper) || input_error("variable $(v.name)", "bounds", "optimise! needs finite bounds")
    end
    free = findall(v -> v.lower < v.upper, vars)
    isempty(free) && input_error("optimise!", "variables", "no variable is free to change")
    lo = [v.lower for v in vars]
    hi = [v.upper for v in vars]
    for (k, v) in enumerate(vars)
        k in free || (p[v.name] = v.lower)
    end
    x0 = clamp.(variable_values(p, names), lo, hi)
    to_x(u) = (x = copy(x0); x[free] .= lo[free] .+ u .* (hi[free] .- lo[free]); x)
    u0 = (x0[free] .- lo[free]) ./ (hi[free] .- lo[free])
    reset_step_pattern!(p)
    sc(k) = Float64(get(scale, k, 1.0))
    history = Any[]
    best = Ref{Any}(nothing)
    function f(u)
        x = to_x(u)
        set_variables!(p, names, x)
        r = performance_gradient(p; variables = names)
        m, g = r.metrics, r.gradients
        grad(k) = [g[k][n] for n in names]
        v = sgn * m[objective] / sc(objective)
        dv = sgn .* grad(objective) ./ sc(objective)
        for (k, target) in pairs(at_least)
            short = (target - m[k]) / sc(k)
            short > 0 && (v -= weight * short^2; dv = dv .+ 2weight * short .* grad(k) ./ sc(k))
        end
        for (k, target) in pairs(at_most)
            over = (m[k] - target) / sc(k)
            over > 0 && (v -= weight * over^2; dv = dv .- 2weight * over .* grad(k) ./ sc(k))
        end
        values = NamedTuple{Tuple(names)}(Tuple(x))
        push!(history, (values = values, metrics = m, objective = v))
        if isnothing(best[]) || v > best[].objective
            best[] = (values = values, metrics = m, objective = v)
        end
        isnothing(log) || log(values, r)
        return (v, dv[free] .* (hi[free] .- lo[free]))
    end
    kw = Dict{Symbol, Any}()
    if !isempty(linear_constraints)
        A = zeros(length(linear_constraints), length(free))
        bvec = zeros(length(linear_constraints))
        for (i, (coeffs, bound)) in enumerate(linear_constraints)
            rhs = Float64(bound)
            for (n, a) in pairs(coeffs)
                j = findfirst(==(Symbol(n)), names)
                isnothing(j) && input_error("optimise!", "linear_constraints", "$n is not among the variables")
                rhs -= a * lo[j]
                jf = findfirst(==(j), free)
                isnothing(jf) || (A[i, jf] = a * (hi[j] - lo[j]))
            end
            bvec[i] = rhs
        end
        kw[:lin_ineq] = (A = A, b = bvec)
    end
    Jutul.unit_box_bfgs(u0, f; maximize = true, max_it = max_it, grad_tol = 1e-4, obj_change_tol = 1e-5,
        max_initial_update = 0.05, stepIncreaseTol = 3.0, print = print_level, kw...)
    set_variables!(p, names, collect(best[].values))
    return (values = best[].values, metrics = best[].metrics, evaluations = length(history), history = history)
end
