"""
    AbstractThermalModel

How a distributed unit treats energy. Choose one when building the system:

- [`WithWall`](@ref): gas-phase energy balance plus a column wall exchanging heat
  with the gas and the surroundings (default)
- [`Adiabatic`](@ref): gas-phase energy balance with no heat exchange through the
  side of the unit
- [`Isothermal`](@ref): temperature held at a fixed value, no energy balance
- [`HeatExchanger`](@ref): gas-phase energy balance plus tubes carrying a heat
  transfer fluid through the unit
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

"""
    HeatExchanger()

Gas-phase energy balance plus tubes running through the unit, carrying a heat
transfer fluid that heats or cools the gas and sorbent indirectly (Mondino et
al. 2022, eqs. A.10–A.11). Each cell gets two more temperatures: the tube
wall, `TubeTemperature`, exchanging heat with the gas and sorbent outside and
with the fluid inside, and the fluid, `HXFluidTemperature`, which flows from
the bottom of the unit to the top.

The fluid enters at the `HXFluidInletTemperature` parameter, or at the outlet
temperature of another unit's tubes when the two are joined by a
`HeatTransferFluidConnection`. Cross-flow tube bundles are treated as
counter-current along the unit, as in the paper.
"""
struct HeatExchanger <: AbstractThermalModel end

"Whether `sys` solves a wall energy balance."
has_wall(sys::MoccaSystem) = thermal_model(sys) isa WithWall

"Whether `sys` solves a gas-phase energy balance (as opposed to a fixed temperature)."
has_energy_balance(sys::MoccaSystem) = !(thermal_model(sys) isa Isothermal)

"Whether `sys` has tubes carrying a heat transfer fluid through it; see [`HeatExchanger`](@ref)."
has_heat_exchanger(sys::MoccaSystem) = thermal_model(sys) isa HeatExchanger

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
