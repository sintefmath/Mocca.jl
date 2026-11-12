"""
    AbstractIsotherm

Abstract base type for adsorption isotherm models.

Implementations must provide:
- `compute_equilibrium(iso, C, T)`: equilibrium loading per component [mol/m³ adsorbent]
- `compute_enthalpy(iso, C, T)`: isosteric heat of adsorption per component [J/mol]
"""
abstract type AbstractIsotherm end

"""
    compute_equilibrium(iso::AbstractIsotherm, C, T) → qstar

Compute equilibrium loading per component [mol/m³ adsorbent] for given concentrations [mol/m³]
and temperature [K]. Returns a vector with one entry per component.
"""
function compute_equilibrium end

"""
    compute_enthalpy(iso::AbstractIsotherm, C, T) → ΔH

Compute isosteric heat of adsorption per component [J/mol] for given concentrations
and temperature.
"""
function compute_enthalpy end

"""
    saturation_capacity(iso::AbstractIsotherm) → capacities or nothing

Saturation capacity per component [mol/m³ adsorbent], the loading the isotherm
tends to at high pressure, or `nothing` if it has none. Optional: it sets the
scale of the loading when changes over a cycle are measured.
"""
saturation_capacity(iso::AbstractIsotherm) = nothing


"""
    minimum_mole_fraction(iso::AbstractIsotherm) → Float64

Smallest mole fraction the gas may take in a bed of this sorbent; by default
10⁻¹². An isotherm that holds a trace component very strongly, so that the
gas in equilibrium with a lightly loaded sorbent holds less of it than that,
needs a smaller one: Newton's method stalls where the fraction it heads for is
below the bound.
"""
minimum_mole_fraction(iso::AbstractIsotherm) = MIN_GAS_MOLEFRACTION
