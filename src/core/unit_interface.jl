# What the flowsheet and process layers need from a unit with holdup.
#
# A unit is a Jutul.SimulationModel of a MoccaSystem. To be part of a
# flowsheet, its system provides `component_names`, and its model provides:
#
# - `ports(model)`: named ports, each the index of the cell a stream enters or
#   leaves through (see coupling/ports.jl);
# - the variables the port cross terms read at a port cell: `Pressure`,
#   `Temperature` and `y`; and, for streams into an energy balance, the
#   system's `heat_capacity_gas` and `molecular_masses`;
# - the `TimeScale` parameter, used through `_scaled_dt` in every time
#   derivative, for gradients with respect to stage durations;
# - [`stream_equations`](@ref): the equations a stream adds to;
# - [`setup_unit_state`](@ref) and [`setup_unit_parameters`](@ref), if the
#   defaults below do not fit.
#
# A DistributedUnit gets the gas block for free (blocks/gas.jl): its face flux
# reads the `Transmissibilities`, `DispersionTransmissibilities`,
# `AxialDispersion` and `FluidViscosity` parameters. A unit that also
# [`holds_sorbent`](@ref) counts as a bed: steps can act on it and its sorbent
# counts in the productivity of a process.

"""
    stream_equations(model) → Tuple of Symbols

Equations of `model` that a stream through one of its ports adds to: the mass
balance, and the energy balance if the unit solves one.
"""
function stream_equations end

function stream_equations(model::DistributedUnitModel)
    return has_energy_balance(model.system) ? (:mass_conservation, :energy_column) : (:mass_conservation,)
end

function stream_equations(model::MoccaModel)
    error("$(nameof(typeof(model.system))) does not define Mocca.stream_equations, so streams cannot connect to it")
end

"""
    setup_unit_parameters(model; overrides...) → Dict

Parameters of a unit, with `overrides` for any of them.
"""
# Typed on MoccaModel rather than Jutul.SimulationModel: Julia ranks a plain
# SimulationModel above SimulationModel{<:Any, <:SorbentBed, <:Any, <:Any}, so
# the method for beds would never be called
setup_unit_parameters(model::MoccaModel; kwarg...) = Jutul.setup_parameters(model; kwarg...)

"""
    setup_unit_state(model; pressure, temperature, composition,
        wall_temperature = temperature) → Dict

Uniform initial state of a unit with holdup. The default sets the gas state,
and the wall temperature if the unit has a wall; a unit with more primary
variables adds a method.
"""
function setup_unit_state(model::MoccaModel; pressure, temperature, composition, wall_temperature = temperature)
    return Jutul.setup_state(model; _gas_state(model, pressure, temperature, composition, wall_temperature)...)
end

function _gas_state(model, pressure, temperature, composition, wall_temperature)
    kw = Dict{Symbol, Any}(:Pressure => pressure, :Temperature => temperature, :y => collect(composition))
    haskey(model.primary_variables, :WallTemperature) && (kw[:WallTemperature] = wall_temperature)
    return kw
end

"""
    holds_sorbent(model) → Bool

Whether `model` is a bed: a unit holding sorbent, which steps act on and whose
sorbent counts in the productivity of a process.
"""
holds_sorbent(model) = false
holds_sorbent(::SorbentBedModel) = true
