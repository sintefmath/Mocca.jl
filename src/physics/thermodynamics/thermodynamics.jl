"""
    AbstractGasThermodynamics

Abstract base type for thermodynamic models of a gas mixture: the heat
capacity and enthalpy of each component, and the equation of state.

Implementations must provide:
- `component_heat_capacities(model, T)`: molar heat capacity of each component
  [J/(mol·K)] at temperature `T` [K], as an `SVector`.

and may provide:
- `component_enthalpies(model, T)`: molar enthalpy of each component [J/mol]
  relative to [`THERMODYNAMIC_REFERENCE_TEMPERATURE`](@ref); by default the
  heat capacities integrated by Gauss–Legendre quadrature;
- `molar_concentration(model, P, T)`: by default that of the ideal gas. The
  energy balances of Mocca assume an ideal gas.

The mixture properties ([`molar_heat_capacity`](@ref), [`molar_enthalpy`](@ref),
[`heat_capacity_ratio`](@ref)) follow from these. Every function must work
for automatic-differentiation numbers and not allocate, as it is evaluated in
every cell at every Newton iteration.
"""
abstract type AbstractGasThermodynamics end

"Reference temperature [K] of [`component_enthalpies`](@ref), at which they are zero."
const THERMODYNAMIC_REFERENCE_TEMPERATURE = 298.15

"""
    component_heat_capacities(model::AbstractGasThermodynamics, T) → SVector

Molar heat capacity [J/(mol·K)] of each component at temperature `T` [K].
"""
function component_heat_capacities end

Framework.number_of_components(model::AbstractGasThermodynamics) =
    length(component_heat_capacities(model, THERMODYNAMIC_REFERENCE_TEMPERATURE))

# 4-point Gauss–Legendre rule on [-1, 1]: exact for heat capacities that are
# polynomials of degree 7 or less
const _GL_NODES = (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const _GL_WEIGHTS = (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)

"""
    component_enthalpies(model::AbstractGasThermodynamics, T) → SVector

Molar enthalpy [J/mol] of each component at temperature `T` [K], relative to
[`THERMODYNAMIC_REFERENCE_TEMPERATURE`](@ref): ``h_i(T) = \\int_{T_{ref}}^T c_{p,i}\\,dT``.
"""
function component_enthalpies(model::AbstractGasThermodynamics, T)
    T0 = THERMODYNAMIC_REFERENCE_TEMPERATURE
    half, mid = (T - T0) / 2, (T + T0) / 2
    cp(k) = _GL_WEIGHTS[k] * component_heat_capacities(model, mid + half * _GL_NODES[k])
    return half * (cp(1) + cp(2) + cp(3) + cp(4))
end

# Σᵢ yᵢ·vᵢ, without allocating
@inline function _mole_average(y, v)
    s = zero(y[1] * v[1])
    for i in eachindex(v)
        s += y[i] * v[i]
    end
    return s
end

"""
    molar_heat_capacity(model::AbstractGasThermodynamics, y, T)

Molar heat capacity [J/(mol·K)] of the mixture of mole fractions `y` at
temperature `T` [K]: ``c_p = \\sum_i y_i c_{p,i}(T)``, the mole average of the
molar heat capacities, or ``\\sum_i y_i M_i c_{p,i}`` with the heat capacities
per mass.
"""
molar_heat_capacity(model::AbstractGasThermodynamics, y, T) = _mole_average(y, component_heat_capacities(model, T))

"""
    molar_enthalpy(model::AbstractGasThermodynamics, y, T)

Molar enthalpy [J/mol] of the mixture of mole fractions `y` at temperature `T`
[K], relative to [`THERMODYNAMIC_REFERENCE_TEMPERATURE`](@ref).
"""
molar_enthalpy(model::AbstractGasThermodynamics, y, T) = _mole_average(y, component_enthalpies(model, T))

"""
    molar_concentration(model::AbstractGasThermodynamics, P, T)

Total molar concentration [mol/m³] of the gas at pressure `P` [Pa] and
temperature `T` [K]: ``P/(RT)`` for the ideal gas.
"""
molar_concentration(model::AbstractGasThermodynamics, P, T) = P / (GAS_CONSTANT * T)

"""
    heat_capacity_ratio(model::AbstractGasThermodynamics, y, T)

Ratio of the heat capacities at constant pressure and volume, ``\\gamma =
c_p/(c_p - R)``, of the ideal gas mixture of mole fractions `y` at temperature
`T` [K].
"""
function heat_capacity_ratio(model::AbstractGasThermodynamics, y, T)
    cp = molar_heat_capacity(model, y, T)
    return cp / (cp - GAS_CONSTANT)
end
