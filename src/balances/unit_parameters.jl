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

# Every distributed unit takes the jacket temperature and the time scale as
# forces. The jacket temperature comes first: Jutul's generic loops look the
# length of a force up by a counter that skips absent forces, and the time
# scale is only set after the forces are devectorised (see force_gradients).
Jutul.setup_forces(model::DistributedUnitModel; jacket_temperature = nothing, time_scale = nothing) =
    (jacket_temperature = jacket_temperature, time_scale = time_scale)

# The time step of a unit's time derivatives
@inline _scaled_dt(state, Δt) = _state_has(state, :TimeScale) ? Δt * state.TimeScale[1] : Δt

"""
    JacketTemperature

Temperature [K] of the medium outside the wall of a unit, such as the heating
or cooling fluid of a jacket, which the wall exchanges heat with through its
outer surface and its two ends. A stage sets it through the unit's
`jacket_temperature` force, which a step's `jacket_temperature` gives; in a
stage that does not, it is the unit's `AmbientTemperature`. Like the time
scale, it is a parameter copied from the forces before each step.
"""
struct JacketTemperature <: Jutul.ScalarVariable end
Jutul.associated_entity(::JacketTemperature) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::JacketTemperature, symb)
    return copy(data_domain[:ambient_temperature, Unit()])
end

function Jutul.update_parameter_before_step!(x, ::JacketTemperature, storage, model, dt, forces)
    T = hasproperty(forces, :jacket_temperature) ? forces.jacket_temperature : nothing
    x .= isnothing(T) ? storage.state.AmbientTemperature : T
    return x
end

# A jacket temperature is a number among the forces of a stage, with a
# gradient. The loops over the forces are those of a device (see
# devices/device.jl), since Jutul's generic ones assume that the unit's values
# start at the beginning of the whole vector of a MultiModel.
Jutul.vectorize_force!(v, model::DistributedUnitModel, f::Real, name, variant) = (v[1] = f; (names = [name],))
Jutul.devectorize_force(f::Real, model::DistributedUnitModel, X, meta, name, variant) = X[1]

function Jutul.vectorize_forces!(v, model::DistributedUnitModel, config, forces; update_config = false)
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

function Jutul.devectorize_forces(forces, model::DistributedUnitModel, X, config; offset = 0)
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
