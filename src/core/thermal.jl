"""
    AbstractThermalModel

How a distributed unit treats energy. Choose one when building the system:

- [`WithWall`](@ref): gas-phase energy balance plus a column wall exchanging heat
  with the gas and the surroundings (default)
- [`Adiabatic`](@ref): gas-phase energy balance with no heat exchange through the
  side of the unit
- [`Isothermal`](@ref): temperature held at a fixed value, no energy balance
"""
abstract type AbstractThermalModel end

"Gas-phase energy balance plus a wall energy balance, with heat exchange between them and to ambient."
struct WithWall <: AbstractThermalModel end

"Gas-phase energy balance without heat exchange through the side of the unit."
struct Adiabatic <: AbstractThermalModel end

"""
    Isothermal(T)

Temperature held at `T` [K] throughout the unit. The value is stored as the
`IsothermalTemperature` parameter, so it can be changed per case and
differentiated.
"""
struct Isothermal{T<:Real} <: AbstractThermalModel
    temperature::T
end

"Whether `sys` solves a wall energy balance."
has_wall(sys::MoccaSystem) = thermal_model(sys) isa WithWall

"Whether `sys` solves a gas-phase energy balance (as opposed to a fixed temperature)."
has_energy_balance(sys::MoccaSystem) = !(thermal_model(sys) isa Isothermal)

"The [`AbstractThermalModel`](@ref) of `sys`."
thermal_model(sys::MoccaSystem) = sys.thermal

# Temperature fixed by a parameter, used in place of the energy balance for
# isothermal units.
struct FixedTemperatureEquation <: Jutul.JutulEquation end

Jutul.associated_entity(::FixedTemperatureEquation) = Jutul.Cells()
Jutul.number_of_equations_per_entity(model::Jutul.SimulationModel, ::FixedTemperatureEquation) = 1
Jutul.local_discretization(::FixedTemperatureEquation, i) = nothing

function Jutul.update_equation_in_entity!(eq_buf, self_cell, state, state0, eq::FixedTemperatureEquation, model, Δt, ldisc = nothing)
    eq_buf[1] = state.Temperature[self_cell] - state.IsothermalTemperature[1]
end

struct IsothermalTemperature <: Jutul.ScalarVariable end
Jutul.associated_entity(::IsothermalTemperature) = Unit()

function Jutul.default_parameter_values(data_domain, model, ::IsothermalTemperature, symb)
    return [thermal_model(model.system).temperature]
end
