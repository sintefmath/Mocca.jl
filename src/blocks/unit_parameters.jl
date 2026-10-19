# Unit-entity Jutul parameters
# Each is a single scalar associated with the Unit entity (count=1),
# accessed at runtime as state.ParamName[1].

"""
    TimeScale

Factor on the time step in every time derivative of a unit. It is 1 in a
simulation. To differentiate with respect to a stage's duration,
[`force_gradients`](@ref) sets it for each step, through the unit's
`time_scale` force, to the stage's relative duration, so that stretching a
stage stretches its steps.
"""
struct TimeScale <: Jutul.ScalarVariable end
Jutul.associated_entity(::TimeScale) = Unit()
Jutul.default_parameter_values(data_domain, model, ::TimeScale, symb) = [1.0]

# A parameter is set once per simulation, while the time scale differs per
# stage, so it comes in as a force and is copied before each step
function Jutul.update_parameter_before_step!(x, ::TimeScale, storage, model, dt, forces)
    s = _time_scale_force(forces)
    isnothing(s) || (x .= s)
    return x
end
_time_scale_force(forces) = hasproperty(forces, :time_scale) ? forces.time_scale : nothing

# The time step of a unit's time derivatives
@inline _scaled_dt(state, Δt) = hasproperty(state, :TimeScale) ? Δt * state.TimeScale[1] : Δt

struct AdsorbentDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::AdsorbentDensity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::AdsorbentDensity, symb)
    return copy(data_domain[:adsorbent_density, Unit()])
end

struct AdsorbentHeatCapacity <: Jutul.ScalarVariable end
Jutul.associated_entity(::AdsorbentHeatCapacity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::AdsorbentHeatCapacity, symb)
    return copy(data_domain[:adsorbent_heat_capacity, Unit()])
end

struct WallDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallDensity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::WallDensity, symb)
    return copy(data_domain[:wall_density, Unit()])
end

struct WallHeatCapacity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallHeatCapacity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::WallHeatCapacity, symb)
    return copy(data_domain[:wall_heat_capacity, Unit()])
end

struct WallConductivity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallConductivity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::WallConductivity, symb)
    return copy(data_domain[:wall_conductivity, Unit()])
end

struct FluidViscosity <: Jutul.ScalarVariable end
Jutul.associated_entity(::FluidViscosity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::FluidViscosity, symb)
    return copy(data_domain[:fluid_viscosity, Unit()])
end

struct FluidDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::FluidDensity) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::FluidDensity, symb)
    return copy(data_domain[:fluid_density, Unit()])
end

struct InnerHeatTransferCoeff <: Jutul.ScalarVariable end
Jutul.associated_entity(::InnerHeatTransferCoeff) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::InnerHeatTransferCoeff, symb)
    return copy(data_domain[:inner_htc, Unit()])
end

struct OuterHeatTransferCoeff <: Jutul.ScalarVariable end
Jutul.associated_entity(::OuterHeatTransferCoeff) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::OuterHeatTransferCoeff, symb)
    return copy(data_domain[:outer_htc, Unit()])
end

struct AmbientTemperature <: Jutul.ScalarVariable end
Jutul.associated_entity(::AmbientTemperature) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::AmbientTemperature, symb)
    return copy(data_domain[:ambient_temperature, Unit()])
end

struct WallCrossSectionArea <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallCrossSectionArea) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::WallCrossSectionArea, symb)
    r_in = first(data_domain[:r_in, Unit()])
    r_out = first(data_domain[:r_out, Unit()])
    return [π * (r_out^2 - r_in^2)]
end
