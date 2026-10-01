# Flow devices: valves, flow controllers and (later) pumps and compressors.
#
# A device is a one-cell model that owns the molar flow F through it, from its
# inlet port to its outlet port, and holds a copy of the pressure, temperature
# and composition at each port. Keeping everything in one cell means a device
# law only depends on the device's own variables, so it can be switched per
# stage through forces while keeping exact derivatives.
#
# Each port is either connected to a unit (see coupling/), which sets the port
# copy through a cross term, or is a boundary whose state is prescribed by a
# PortCondition force.

"""
    FlowDevice(component_names)

A flow element with an inlet and an outlet port. The flow law (closed, valve,
flow controller, ...) is set per stage through the `law` force; see
[`AbstractDeviceLaw`](@ref). Build its model with [`setup_device_model`](@ref).
"""
struct FlowDevice <: FlowElement
    component_names::Vector{String}
end

const FlowDeviceModel = Jutul.SimulationModel{<:Jutul.JutulDomain, <:FlowDevice}

# Residual scales: the port equations are linear, the law is mildly nonlinear.
# Scaling the flow residual to mmol/s and fractions by 1e3 makes the default
# nonlinear tolerance (1e-3) tight for these equations.
const DEVICE_FLOW_SCALE = 1e3
const DEVICE_FRACTION_SCALE = 1e3

"Molar flow [mol/s] through a device, positive from inlet to outlet."
struct DeviceMolarFlow <: Jutul.ScalarVariable end
Jutul.minimum_value(::DeviceMolarFlow) = -Inf

"""
    setup_device_model(component_names)

Build the `SimulationModel` for a [`FlowDevice`](@ref).
"""
function setup_device_model(component_names)
    mesh = Jutul.CartesianMesh((1, 1, 1), (1.0, 1.0, 1.0))
    domain = Jutul.DataDomain(mesh)
    domain.entities[Unit()] = 1
    return Jutul.SimulationModel(domain, FlowDevice(collect(String, component_names)))
end

"Port names of a device, in order: inlet (1) and outlet (2)."
const DEVICE_PORTS = (:inlet, :outlet)

port_index(::FlowDeviceModel, port::Symbol) = port == :inlet ? 1 : (port == :outlet ? 2 : error("Unknown device port $port, expected :inlet or :outlet"))

const _PORT_VARIABLES = (
    (pressure = :InletPressure, temperature = :InletTemperature, composition = :InletComposition),
    (pressure = :OutletPressure, temperature = :OutletTemperature, composition = :OutletComposition),
)

@inline port_pressure(state, k) = k == 1 ? state.InletPressure[1] : state.OutletPressure[1]
@inline port_temperature(state, k) = k == 1 ? state.InletTemperature[1] : state.OutletTemperature[1]
@inline port_composition(state, k, i) = k == 1 ? state.InletComposition[i, 1] : state.OutletComposition[i, 1]

function Jutul.select_primary_variables!(S, system::FlowDevice, model::Jutul.SimulationModel)
    # Port copies are set by linear equations, so they need no increment limits
    for pv in _PORT_VARIABLES
        S[pv.pressure] = Pressure(max_rel = Inf)
        S[pv.temperature] = Temperature(max_rel = Inf)
        S[pv.composition] = GasMoleFractions(dz_max = Inf)
    end
    S[:MolarFlow] = DeviceMolarFlow()
end

Jutul.select_secondary_variables!(S, system::FlowDevice, model::Jutul.SimulationModel) = nothing
Jutul.select_parameters!(S, system::FlowDevice, model::Jutul.SimulationModel) = nothing

# Port state equation: the port copy equals whatever it is attached to,
#   [P_k, T_k, s·y_k,1 … s·y_k,N−1] − [P, T, s·y_1 … s·y_N−1]_attached = 0
# The local part is the port copy; a PortStateCT or a boundary force subtracts
# the attached state.
struct PortStateEquation{K} <: Jutul.JutulEquation end

Jutul.associated_entity(::PortStateEquation) = Jutul.Cells()
Jutul.local_discretization(::PortStateEquation, i) = nothing
Jutul.number_of_equations_per_entity(model::FlowDeviceModel, ::PortStateEquation) = number_of_components(model.system) + 1

function Jutul.update_equation_in_entity!(eq_buf, self_cell, state, state0, eq::PortStateEquation{K}, model, Δt, ldisc = nothing) where K
    N = number_of_components(model.system)
    eq_buf[1] = port_pressure(state, K)
    eq_buf[2] = port_temperature(state, K)
    for i in 1:(N - 1)
        eq_buf[2 + i] = DEVICE_FRACTION_SCALE * port_composition(state, K, i)
    end
end

# Device law equation: s·(F − F_law) = 0. The local part is s·F; the law force
# subtracts s·F_law. Without a law force the device is closed.
struct DeviceLawEquation <: Jutul.JutulEquation end

Jutul.associated_entity(::DeviceLawEquation) = Jutul.Cells()
Jutul.local_discretization(::DeviceLawEquation, i) = nothing
Jutul.number_of_equations_per_entity(model::FlowDeviceModel, ::DeviceLawEquation) = 1

function Jutul.update_equation_in_entity!(eq_buf, self_cell, state, state0, eq::DeviceLawEquation, model, Δt, ldisc = nothing)
    eq_buf[1] = DEVICE_FLOW_SCALE * state.MolarFlow[1]
end

function Jutul.select_equations!(eqs, system::FlowDevice, model::Jutul.SimulationModel)
    eqs[:inlet_state] = PortStateEquation{1}()
    eqs[:outlet_state] = PortStateEquation{2}()
    eqs[:device_law] = DeviceLawEquation()
end

"Equation holding the state of device port `port` (`:inlet` or `:outlet`)."
port_equation(port::Symbol) = port == :inlet ? :inlet_state : :outlet_state

# ------------------------------------------------------------------------------
# Device laws
# ------------------------------------------------------------------------------

"""
    AbstractDeviceLaw

How a [`FlowDevice`](@ref) relates its flow to the states at its ports. Passed
as the `law` force for each stage. Implementations provide
[`device_flow`](@ref).
"""
abstract type AbstractDeviceLaw end

"""
    device_flow(law, state, time) → F

Molar flow [mol/s] from inlet to outlet prescribed by `law`, given the device
`state` (with port pressures, temperatures and compositions).
"""
function device_flow end

"No flow."
struct Closed <: AbstractDeviceLaw end
device_flow(::Closed, state, time) = false

"""
    LinearValve(conductance; opening = 1.0)

Valve whose volumetric flow at upstream conditions is linear in the pressure
drop, `q = opening·conductance·(P_in − P_out)`, with conductance in
m³/(s·Pa). The molar flow uses the upstream concentration, so it changes sign
smoothly through zero.
"""
struct LinearValve{T<:Real} <: AbstractDeviceLaw
    conductance::T
    opening::T
end
LinearValve(conductance; opening = 1.0) = LinearValve(promote(conductance, opening)...)

function device_flow(law::LinearValve, state, time)
    P_in, P_out = port_pressure(state, 1), port_pressure(state, 2)
    ΔP = P_in - P_out
    up = ΔP >= 0 ? 1 : 2
    c_up = port_pressure(state, up) / (GAS_CONSTANT * port_temperature(state, up))
    return law.opening * law.conductance * ΔP * c_up
end

"""
    VolumetricFlow(rate)

Flow controller delivering `rate` [m³/s] of gas, measured at the outlet
pressure and inlet temperature, from inlet to outlet. Use it for a feed
specified as a superficial velocity at the bed inlet: `rate = v·A`.
"""
struct VolumetricFlow{T<:Real} <: AbstractDeviceLaw
    rate::T
end

function device_flow(law::VolumetricFlow, state, time)
    return law.rate * port_pressure(state, 2) / (GAS_CONSTANT * port_temperature(state, 1))
end

function Jutul.apply_forces_to_equation!(acc, storage, model::FlowDeviceModel, eq::DeviceLawEquation, eq_s, law::AbstractDeviceLaw, time)
    isnothing(acc) && return  # parameter sensitivities, see below
    acc[1] -= DEVICE_FLOW_SCALE * device_flow(law, storage.state, time)
end

# When Jutul computes sensitivities with respect to parameters it swaps each
# model's primary variables for its parameters. A device has no parameters, so
# its equations have no diagonal entries then and forces get `nothing`; the
# device equations do not depend on any parameter, so there is nothing to add.

# ------------------------------------------------------------------------------
# Port conditions (boundaries)
# ------------------------------------------------------------------------------

"""
    PortCondition(; pressure, temperature, composition)

Prescribed state at a device port that is not connected to a unit, i.e. a
boundary such as a feed, a product header or a vacuum line. `pressure` may be a
number or a time-dependent profile such as [`ExponentialRamp`](@ref).
"""
struct PortCondition{P, T<:Real, Y<:AbstractVector}
    pressure::P
    temperature::T
    composition::Y
end
PortCondition(; pressure, temperature, composition) = PortCondition(pressure, temperature, composition)

"""
    ExponentialRamp(start, stop, λ; t0 = 0.0)

Value moving from `start` towards `stop` as `stop + (start − stop)·exp(−λ(t − t0))`.
The pressure profile used for pressurisation, blowdown and evacuation stages.
"""
struct ExponentialRamp{T<:Real}
    start::T
    stop::T
    λ::T
    t0::Float64
end
ExponentialRamp(start, stop, λ; t0 = 0.0) = ExponentialRamp(promote(start, stop, λ)..., Float64(t0))

evaluate_profile(x::Real, time) = x
evaluate_profile(r::ExponentialRamp, time) = r.stop + (r.start - r.stop) * exp(-r.λ * (time - r.t0))

"Shift a time-dependent profile so that its local time starts at `t0`."
shift_profile(x::Real, t0) = x
shift_profile(r::ExponentialRamp, t0) = ExponentialRamp(r.start, r.stop, r.λ, r.t0 + t0)
shift_profile(c::PortCondition, t0) = PortCondition(shift_profile(c.pressure, t0), c.temperature, c.composition)

# Boundary force on a specific port
struct PortSetpoint{K, C<:PortCondition}
    condition::C
end
PortSetpoint{K}(c::PortCondition) where K = PortSetpoint{K, typeof(c)}(c)

function Jutul.apply_forces_to_equation!(acc, storage, model::FlowDeviceModel, eq::PortStateEquation{K}, eq_s, force::PortSetpoint{K}, time) where K
    isnothing(acc) && return  # parameter sensitivities, see the device law force
    c = force.condition
    N = number_of_components(model.system)
    acc[1] -= evaluate_profile(c.pressure, time)
    acc[2] -= c.temperature
    for i in 1:(N - 1)
        acc[2 + i] -= DEVICE_FRACTION_SCALE * c.composition[i]
    end
end

"""
    Jutul.setup_forces(model::FlowDeviceModel; law = Closed(), inlet = nothing, outlet = nothing)

Forces for one stage of a device: the flow `law`, and a [`PortCondition`](@ref)
for each port that is a boundary.
"""
function Jutul.setup_forces(model::FlowDeviceModel; law::AbstractDeviceLaw = Closed(), inlet = nothing, outlet = nothing)
    wrap(c, k) = (isnothing(c) || c isa PortSetpoint) ? c : PortSetpoint{k}(c)
    return (law = law, inlet = wrap(inlet, 1), outlet = wrap(outlet, 2))
end

# ------------------------------------------------------------------------------
# Force vectorisation, for gradients with respect to device settings
# ------------------------------------------------------------------------------
#
# Jutul.solve_adjoint_forces differentiates with respect to the numbers that
# make up the forces of each stage:
#
#   law      LinearValve: conductance, opening; VolumetricFlow: rate; Closed: none
#   inlet,   pressure (a number, or start, stop and λ of an ExponentialRamp),
#   outlet   temperature, and each mole fraction
#
# Mole fractions are treated as independent, so their gradients do not account
# for the fractions summing to one.

_profile_values(x::Real) = [x]
_profile_values(r::ExponentialRamp) = [r.start, r.stop, r.λ]
_profile_names(::Real) = [:pressure]
_profile_names(::ExponentialRamp) = [:pressure_start, :pressure_stop, :pressure_λ]
_profile_from(::Real, X) = X[1]
_profile_from(r::ExponentialRamp, X) = ExponentialRamp(X[1], X[2], X[3]; t0 = r.t0)

_force_values(::Closed) = Float64[]
_force_values(law::LinearValve) = [law.conductance, law.opening]
_force_values(law::VolumetricFlow) = [law.rate]
_force_values(f::PortSetpoint) = vcat(_profile_values(f.condition.pressure), f.condition.temperature, f.condition.composition)

_force_names(::Closed) = Symbol[]
_force_names(::LinearValve) = [:conductance, :opening]
_force_names(::VolumetricFlow) = [:rate]
_force_names(f::PortSetpoint) = vcat(_profile_names(f.condition.pressure), :temperature,
    [Symbol("y_$i") for i in eachindex(f.condition.composition)])

_force_from(::Closed, X) = Closed()
_force_from(::LinearValve, X) = LinearValve(X[1], X[2])
_force_from(::VolumetricFlow, X) = VolumetricFlow(X[1])
function _force_from(f::PortSetpoint{K}, X) where K
    c = f.condition
    n = length(_profile_values(c.pressure))
    return PortSetpoint{K}(PortCondition(_profile_from(c.pressure, X[1:n]), X[n + 1], collect(X[(n + 2):end])))
end

const _DeviceForce = Union{AbstractDeviceLaw, PortSetpoint}

Jutul.vectorization_length(f::_DeviceForce, model::FlowDeviceModel, name, variant) = length(_force_values(f))

function Jutul.vectorize_force!(v, model::FlowDeviceModel, f::_DeviceForce, name, variant)
    v .= _force_values(f)
    return (names = _force_names(f),)
end

Jutul.devectorize_force(f::_DeviceForce, model::FlowDeviceModel, X, meta, name, variant) = _force_from(f, X)

# Jutul's generic loops over a model's forces store a length for every force
# but look the lengths up by a counter that skips absent forces, and a
# MultiModel passes each submodel the whole vector. A device's law is always
# present while its port conditions often are not, so devices do the loops
# themselves, starting at the device's offset within the vector.
function Jutul.vectorize_forces!(v, model::FlowDeviceModel, config, forces; update_config = false)
    offset = first(config.offsets) - 1
    for (i, (k, f)) in enumerate(pairs(forces))
        n = config.lengths[i]
        (isnothing(f) || isnothing(config.targets[k])) && continue
        m = Jutul.vectorize_force!(view(v, offset .+ (1:n)), model, f, k, config.targets[k])
        update_config && (config.meta[k] = m)
        offset += n
    end
    return v
end

function Jutul.devectorize_forces(forces, model::FlowDeviceModel, X, config; offset = 0)
    offset = 0
    out = Dict{Symbol, Any}()
    for (i, (k, f)) in enumerate(pairs(forces))
        n = config.lengths[i]
        if isnothing(f) || isnothing(config.targets[k])
            out[k] = f
        else
            out[k] = Jutul.devectorize_force(f, model, view(X, offset .+ (1:n)), config.meta[k], k, config.targets[k])
            offset += n
        end
    end
    return Jutul.setup_forces(model; out...)
end

"""
    setup_device_state(model; inlet, outlet, flow = 0.0)

Initial state of a device from the [`PortCondition`](@ref) (or any object with
`pressure`, `temperature` and `composition`) at each port.
"""
function setup_device_state(model::FlowDeviceModel; inlet, outlet, flow = 0.0)
    return Jutul.setup_state(model;
        InletPressure = evaluate_profile(inlet.pressure, 0.0),
        InletTemperature = inlet.temperature,
        InletComposition = collect(inlet.composition),
        OutletPressure = evaluate_profile(outlet.pressure, 0.0),
        OutletTemperature = outlet.temperature,
        OutletComposition = collect(outlet.composition),
        MolarFlow = flow,
    )
end

Jutul.convergence_criterion(model::FlowDeviceModel, storage, eq::PortStateEquation, eq_s, r; dt = 1, update_report = missing) =
    (AbsMax = (errors = vec(maximum(abs, r; dims = 2)), names = ["P", "T", map(i -> "y_$i", 1:(size(r, 1) - 2))...]),)
