"""
    AbstractIsotherm

Abstract base type for adsorption isotherm models.

Implementations must provide:
- `compute_equilibrium(iso, C, T)`: equilibrium loading per component [mol/kg]
- `compute_enthalpy(iso, C, T)`: isosteric heat of adsorption per component [J/mol]
"""
abstract type AbstractIsotherm end

"""
    compute_equilibrium(iso::AbstractIsotherm, C, T) → qstar

Compute equilibrium loading per component [mol/kg] for given concentrations [mol/m³]
and temperature [K]. Returns a vector with one entry per component.
"""
function compute_equilibrium end

"""
    compute_enthalpy(iso::AbstractIsotherm, C, T) → ΔH

Compute isosteric heat of adsorption per component [J/mol] for given concentrations
and temperature.
"""
function compute_enthalpy end

include("dual_site_langmuir.jl")
