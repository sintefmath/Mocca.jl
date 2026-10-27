# Performance measures as functions of integrals.
#
# Every measure of a process's performance is a function of a few integrals
# over its cycle, such as the moles of the key component in the extract or the
# work of the pumps, and of quantities that do not depend on the state, such as
# the cycle time or the adsorbent volume. At steady state the integrals are
# taken over one second, so they are rates. Each integral is a Jutul sum
# objective, so the adjoints give its gradient, and the chain rule from the
# integrals to a measure follows by forward-mode automatic differentiation of
# the measure (measure_partials).

"""
    SumObjective(G...)

Objective for adjoint gradients: the sum of the objectives `G`, as one
objective. With none, it is zero.
"""
struct SumObjective{T}
    terms::T
end
SumObjective(G...) = SumObjective(G)

(G::SumObjective)(model, state, dt, step_info, forces) =
    sum((g(model, state, dt, step_info, forces) for g in G.terms); init = 0.0)

"""
    BlowerObjective(device; efficiency = 0.72, γ = 1.4, scale = 1.0)

Objective for adjoint gradients: `scale` times the work [J] of a fan or blower
that moves the gas through `device`, a link between two units, up in pressure:
the gas is compressed adiabatically with isentropic `efficiency` from the
pressure and temperature of the upstream port to the pressure of the
downstream one, whichever way it flows. A link that lets the gas down in
pressure, such as a valve, takes no work.
"""
struct BlowerObjective
    device::Symbol
    efficiency::Float64
    γ::Float64
    scale::Float64
end
BlowerObjective(device::Symbol; efficiency = 0.72, γ = 1.4, scale = 1.0) = BlowerObjective(device, efficiency, γ, scale)

function (G::BlowerObjective)(model, state, dt, step_info, forces)
    s = state[G.device]
    F = s[:MolarFlow][1]
    P_in, P_out = s[:InletPressure][1], s[:OutletPressure][1]
    T_in, T_out = s[:InletTemperature][1], s[:OutletTemperature][1]
    # Both ports are read whichever way the gas flows, so that the objective
    # reads the same variables
    forward = F >= 0
    n = ifelse(forward, F, -F)
    P_up, P_down = ifelse(forward, P_in, P_out), ifelse(forward, P_out, P_in)
    T_up = ifelse(forward, T_in, T_out)
    k = (G.γ - 1) / G.γ
    W = n * GAS_CONSTANT * T_up / k * ((P_down / P_up)^k - 1) / G.efficiency
    return G.scale * dt * ifelse(P_down > P_up, W, zero(W))
end

"""
    blower_energy(states, timesteps, device; efficiency = 0.72, γ = 1.4,
        window = eachindex(states))

Work [J] of a fan or blower moving the gas through link `device` up in
pressure over the steps in `window`, as in [`BlowerObjective`](@ref).
"""
function blower_energy(states, timesteps, device::Symbol; efficiency = 0.72, γ = 1.4, window = eachindex(states))
    G = BlowerObjective(device; efficiency = efficiency, γ = γ)
    return sum((G(nothing, states[k], timesteps[k], nothing, nothing) for k in window); init = 0.0)
end

# The heat transfer fluid of a unit with a heat exchanger: its heat capacity
# rate ṁ·c_p and its inlet and outlet temperatures. The fluid enters at the
# outlet temperature of the unit upstream, if it comes from one, or at the
# unit's HXFluidInletTemperature. The state holds the unit's parameters.
function _hx_fluid(state, unit, upstream)
    s = state[unit]
    mc = s[:HXFluidMassFlow][1] * s[:HXFluidHeatCapacity][1]
    T_in = isnothing(upstream) ? s[:HXFluidInletTemperature][1] : state[upstream][:HXFluidTemperature][end]
    return (mc, T_in, s[:HXFluidTemperature][end])
end

# The unit whose heat transfer fluid feeds that of `unit`, if any
function _hx_upstream(fs::AbstractFlowsheet, unit::Symbol)
    for c in connections(fs)
        c isa HeatTransferFluidConnection && c.downstream == unit && return c.upstream
    end
    return nothing
end

"""
    HeatDutyObjective(unit; upstream = nothing, scale = 1.0)
    HeatDutyObjective(fs, unit; scale = 1.0)

Objective for adjoint gradients: `scale` times the heat [J] that the heat
transfer fluid gives the gas and sorbent of `unit`, `ṁ·c_p·(T_in − T_out)` of
the fluid, as in [`heat_duty`](@ref). The fluid enters at the outlet
temperature of the unit `upstream`, or at the unit's
`HXFluidInletTemperature`; the second form finds `upstream` from the
[`HeatTransferFluidConnection`](@ref)s of flowsheet `fs`.
"""
struct HeatDutyObjective
    unit::Symbol
    upstream::Union{Nothing, Symbol}
    scale::Float64
end
HeatDutyObjective(unit::Symbol; upstream = nothing, scale = 1.0) = HeatDutyObjective(unit, upstream, scale)
HeatDutyObjective(fs::AbstractFlowsheet, unit::Symbol; scale = 1.0) = HeatDutyObjective(unit, _hx_upstream(fs, unit), scale)

function (G::HeatDutyObjective)(model, state, dt, step_info, forces)
    mc, T_in, T_out = _hx_fluid(state, G.unit, G.upstream)
    return G.scale * dt * mc * (T_in - T_out)
end

"""
    EquivalentWorkObjective(unit; ambient_temperature = 298.15, upstream = nothing, scale = 1.0)
    EquivalentWorkObjective(fs, unit; ambient_temperature = 298.15, scale = 1.0)

Objective for adjoint gradients: `scale` times the work [J] equivalent to the
heat that the heat transfer fluid gives `unit`, as in
[`HeatDutyObjective`](@ref): the work a reversible engine could have made
from that heat against `ambient_temperature` `T₀`. Each bit of heat `δQ`
leaves the fluid at its temperature `T` along the way from `T_in` to `T_out`,
and is worth the Carnot factor `(1 − T₀/T) δQ`, which sums to

    W = ṁ·c_p·[(T_in − T_out) − T₀·ln(T_in/T_out)],

the exergy the fluid gives up. Heat at a high temperature is worth more than
the same heat at a low one.
"""
struct EquivalentWorkObjective
    unit::Symbol
    upstream::Union{Nothing, Symbol}
    ambient_temperature::Float64
    scale::Float64
end
EquivalentWorkObjective(unit::Symbol; ambient_temperature = 298.15, upstream = nothing, scale = 1.0) =
    EquivalentWorkObjective(unit, upstream, ambient_temperature, scale)
EquivalentWorkObjective(fs::AbstractFlowsheet, unit::Symbol; ambient_temperature = 298.15, scale = 1.0) =
    EquivalentWorkObjective(unit, _hx_upstream(fs, unit), ambient_temperature, scale)

function (G::EquivalentWorkObjective)(model, state, dt, step_info, forces)
    mc, T_in, T_out = _hx_fluid(state, G.unit, G.upstream)
    return G.scale * dt * mc * ((T_in - T_out) - G.ambient_temperature * log(T_in / T_out))
end

"""
    equivalent_work(fs, state, parameters, unit; ambient_temperature = 298.15) → W

Work [W] equivalent to the heat that the heat transfer fluid gives `unit` in
flowsheet `state`, as in [`EquivalentWorkObjective`](@ref).
"""
function equivalent_work(fs::AbstractFlowsheet, state, parameters, unit::Symbol; ambient_temperature = 298.15)
    has_heat_exchanger(fs[unit].system) || error("Unit $unit has no heat exchanger")
    G = EquivalentWorkObjective(fs, unit; ambient_temperature = ambient_temperature)
    return G(nothing, objective_state(state, parameters), 1.0, nothing, nothing)
end

# ------------------------------------------------------------------------------
# Integrals and measures
# ------------------------------------------------------------------------------

"""
    objective_state(state, parameters)

A flowsheet `state` with the `parameters` of each unit added where it does not
have them, as Jutul gives an objective in its adjoints, so that an objective
can read parameters such as the heat transfer fluid's mass flow.
"""
function objective_state(state, parameters)
    # Jutul's local states for automatic differentiation hold the parameters
    state isa Union{AbstractDict, Jutul.JutulStorage} || return state
    out = Dict{Symbol, Any}()
    for (k, v) in pairs(state)
        p = get(parameters, k, nothing)
        out[k] = p isa AbstractDict && !all(n -> haskey(v, n), keys(p)) ? _with_unit_parameters(v, p) : v
    end
    return out
end

function _with_unit_parameters(v, p)
    d = Dict{Symbol, Any}(p)
    for (n, x) in pairs(v)
        d[n] = x
    end
    return d
end

# An objective given the parameters, for Jutul's generic adjoints, which give
# an objective the state without them. They are fixed there, so they are added
# as values.
struct _WithParameters{G, P}
    G::G
    parameters::P
end
(o::_WithParameters)(model, state, dt, step_info, forces) = o.G(model, objective_state(state, o.parameters), dt, step_info, forces)

"""
    integral_values(integrals, model, states, timesteps, parameters; window = eachindex(states)) → NamedTuple

Value of each objective in the `NamedTuple` `integrals` summed over the steps
in `window`, with the `parameters` added to each state (see
[`objective_state`](@ref)). The objectives get no forces.
"""
function integral_values(integrals::NamedTuple, model, states, timesteps, parameters; window = eachindex(states))
    N = length(window)
    t = cumsum(timesteps[window])
    s = [objective_state(states[k], parameters) for k in window]
    info = [Jutul.optimization_step_info(i, t[i], timesteps[k]; Nstep = N) for (i, k) in enumerate(window)]
    return map(G -> sum((G(model, s[i], timesteps[k], info[i], nothing) for (i, k) in enumerate(window)); init = 0.0), integrals)
end

"""
    measure_partials(f, sums::NamedTuple) → NamedTuple

Partial derivatives of the measure `f(sums)` with respect to each entry of
`sums`, the integrals and quantities it is a function of, by forward-mode
automatic differentiation. The chain rule then gives the gradient of the
measure from those of the integrals.
"""
function measure_partials(f, sums::NamedTuple)
    ks = keys(sums)
    x = collect(Float64, values(sums))
    g = Jutul.ForwardDiff.gradient(v -> f(NamedTuple{ks}(Tuple(v))), x)
    return NamedTuple{ks}(Tuple(g))
end
