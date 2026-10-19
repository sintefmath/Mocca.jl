"""
    restart_state(model, state) → Dict

The primary variables of `state` (for example the last state of a simulation),
set up as a state that can start a new simulation of `model`. Works for a
single unit and for a flowsheet `MultiModel`.
"""
function restart_state(model::Jutul.SimulationModel, state)
    primary = Dict{Symbol, Any}(k => copy(state[k]) for k in keys(model.primary_variables))
    return Jutul.setup_state(model, primary)
end

function restart_state(model::Jutul.MultiModel, state)
    return Dict{Symbol, Any}(k => restart_state(model.models[k], state[k]) for k in Jutul.submodels_symbols(model))
end

"""
    cycle_change(model, state_start, state_end) → Float64

Largest relative change of any primary variable of a unit with holdup between
the start and the end of a cycle: `max |x_end − x_start| / max |x_start|`, taken
per variable. Devices are skipped, since they hold no inventory.
"""
function cycle_change(model::Jutul.SimulationModel, s0, s1)
    err = 0.0
    for k in keys(model.primary_variables)
        x0, x1 = s0[k], s1[k]
        scale = max(maximum(abs, x0), eps())
        err = max(err, maximum(abs, x1 .- x0) / scale)
    end
    return err
end

cycle_change(model::FlowDeviceModel, s0, s1) = 0.0

function cycle_change(model::Jutul.MultiModel, s0, s1)
    return maximum(cycle_change(model.models[k], s0[k], s1[k]) for k in Jutul.submodels_symbols(model))
end

# ------------------------------------------------------------------------------
# Stacking the state of the units with holdup into one vector, for acceleration
# ------------------------------------------------------------------------------

_holdup_units(model::Jutul.SimulationModel) = [(nothing, model)]
_holdup_units(model::Jutul.MultiModel) = [(k, model.models[k]) for k in Jutul.submodels_symbols(model) if !(model.models[k] isa FlowDeviceModel)]

_unit_state(state, ::Nothing) = state
_unit_state(state, k::Symbol) = state[k]

function _stack(model, state)
    v = Float64[]
    for (k, m) in _holdup_units(model), p in keys(m.primary_variables)
        append!(v, vec(_unit_state(state, k)[p]))
    end
    return v
end

# Scale of each entry: the largest magnitude of its variable in `state`
function _stack_scale(model, state)
    v = Float64[]
    for (k, m) in _holdup_units(model), p in keys(m.primary_variables)
        x = _unit_state(state, k)[p]
        append!(v, fill(max(maximum(abs, x), eps()), length(x)))
    end
    return v
end

_bound(x::Nothing, default) = default
_bound(x::Real, default) = x

# Write `v` into a copy of `state`, keeping each variable within its bounds
function _unstack(model, state, v)
    out = deepcopy(state)
    offset = 0
    for (k, m) in _holdup_units(model), p in keys(m.primary_variables)
        x = _unit_state(out, k)[p]
        n = length(x)
        x .= reshape(v[offset .+ (1:n)], size(x))
        offset += n
        var = m.primary_variables[p]
        x .= clamp.(x, _bound(Jutul.minimum_value(var), -Inf), _bound(Jutul.maximum_value(var), Inf))
        if var isa Jutul.FractionVariables
            x ./= sum(x, dims = 1)
        end
    end
    return out
end

# ------------------------------------------------------------------------------
# Cyclic steady state
# ------------------------------------------------------------------------------

"""
    simulate_to_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-4, max_cycles = 100, acceleration = :anderson, depth = 5,
        info_level = -1, kwarg...)

Repeat one cycle, given by its `forces` and `timesteps`, until the
[`cycle_change`](@ref) of a cycle falls below `tol` or `max_cycles` is reached.
Works for a single bed with legacy boundary conditions and for flowsheets.

Each cycle maps its initial state `x` to its final state `Φ(x)`, and cyclic
steady state is the fixed point `x = Φ(x)`. With `acceleration = :none` the next
cycle starts from `Φ(x)` (successive substitution). With `:anderson` it starts
from a combination of the last `depth` cycles that extrapolates towards the
fixed point (Anderson acceleration), which usually needs far fewer cycles. The
extrapolated state is kept within the bounds of each variable, and if a cycle
started from it fails, the next cycle starts from `Φ(x)` instead.

Returns a NamedTuple with the states and time steps of the last cycle, its
initial state `state0`, the change per cycle in `history`, the number of
`cycles` simulated, and whether it `converged`. Extra keywords go to
[`simulate_process`](@ref).
"""
function simulate_to_cyclic_steady_state(model, state0, parameters, forces, timesteps;
        tol = 1e-4, max_cycles = 100, acceleration = :anderson, depth = 5, info_level = -1, kwarg...)
    acceleration in (:none, :anderson) || error("acceleration must be :none or :anderson, got $acceleration")
    cycle_time = sum(timesteps)
    function run_cycle(x)
        case = MoccaCase(model, timesteps, forces; state0 = x, parameters = parameters)
        states, ts = simulate_process(case; info_level = info_level, output_substates = true, kwarg...)
        ok = !isempty(ts) && isapprox(sum(ts), cycle_time)
        return (states, ts, ok)
    end

    x = restart_state(model, state0)
    history = Float64[]
    # Anderson history, in scaled variables
    ΔF = Vector{Float64}[]
    ΔG = Vector{Float64}[]
    f_prev = g_prev = nothing
    scale = nothing
    # End of the last successful cycle, the fallback start
    last_end = nothing
    states = ts = nothing
    for cycle in 1:max_cycles
        states, ts, ok = run_cycle(x)
        if !ok
            isnothing(last_end) && error("Simulation of cycle $cycle failed")
            empty!(ΔF); empty!(ΔG)
            f_prev = g_prev = nothing
            x = last_end
            states, ts, ok = run_cycle(x)
            ok || error("Simulation of cycle $cycle failed, also from the unaccelerated start")
        end
        x_end = restart_state(model, states[end])
        push!(history, cycle_change(model, x, x_end))
        if history[end] < tol
            return (states = states, timesteps = ts, state0 = x, history = history, cycles = cycle, converged = true)
        end
        last_end = x_end
        if acceleration == :none
            x = x_end
            continue
        end

        isnothing(scale) && (scale = _stack_scale(model, x_end))
        g = _stack(model, x_end) ./ scale
        f = g .- _stack(model, x) ./ scale
        if !isnothing(f_prev)
            push!(ΔF, f .- f_prev)
            push!(ΔG, g .- g_prev)
            if length(ΔF) > depth
                popfirst!(ΔF)
                popfirst!(ΔG)
            end
        end
        f_prev, g_prev = f, g
        if isempty(ΔF)
            x_next = g
        else
            γ = reduce(hcat, ΔF) \ f
            x_next = g .- reduce(hcat, ΔG) * γ
        end
        x = restart_state(model, _unstack(model, x_end, x_next .* scale))
    end
    return (states = states, timesteps = ts, state0 = x, history = history, cycles = max_cycles, converged = false)
end

"""
    simulate_until_mass_balance(model, state0, parameters, forces, timesteps;
        inflow, outflow, component = 1, tol = 0.005, consecutive = 5,
        min_cycles = 50, max_cycles = 300, info_level = -1, kwarg...)

Repeat one cycle of a flowsheet, given by its `forces` and `timesteps`, until
the mass balance of `component` over the flowsheet closes: the relative
difference between the moles entering through the `inflow` devices and leaving
through the `outflow` devices is below `tol` for `consecutive` cycles in a row,
after at least `min_cycles`. This is the cyclic steady state criterion of
Haghpanah et al. (2013) and Krishnamurthy et al. (2014), with CO2 as the
component.

It only checks one component. A more strongly adsorbed component, such as
water, can still be accumulating when it is met; use
[`simulate_to_cyclic_steady_state`](@ref) to require every variable to settle.

Returns a NamedTuple with the states and time steps of the last cycle, its
initial state `state0`, the balance error of each cycle in `history`, the number
of `cycles` simulated, and whether the criterion was met (`converged`). Extra
keywords go to [`simulate_process`](@ref).
"""
function simulate_until_mass_balance(model::Jutul.MultiModel, state0, parameters, forces, timesteps;
        inflow, outflow, component = 1, tol = 0.005, consecutive = 5,
        min_cycles = 50, max_cycles = 300, info_level = -1, kwarg...)
    cycle_time = sum(timesteps)
    x = restart_state(model, state0)
    history = Float64[]
    states = ts = nothing
    x_start = x
    for cycle in 1:max_cycles
        x_start = x
        case = MoccaCase(model, timesteps, forces; state0 = x, parameters = parameters)
        states, ts = simulate_process(case; info_level = info_level, output_substates = true, kwarg...)
        (!isempty(ts) && isapprox(sum(ts), cycle_time)) || error("Simulation of cycle $cycle failed")
        n_in = sum(d -> stream_totals(states, ts, d)[component], inflow)
        n_out = sum(d -> stream_totals(states, ts, d)[component], outflow)
        push!(history, abs(n_in - n_out) / abs(n_in))
        if cycle >= min_cycles && length(history) >= consecutive && all(<(tol), history[end-consecutive+1:end])
            return (states = states, timesteps = ts, state0 = x_start, history = history, cycles = cycle, converged = true)
        end
        x = restart_state(model, states[end])
    end
    return (states = states, timesteps = ts, state0 = x_start, history = history, cycles = max_cycles, converged = false)
end
