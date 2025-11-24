# Abstract isotherm interface for Mocca.jl
# This module provides a modular way to support different isotherm models

"""
    AbstractIsothermModel

Abstract base type for all isotherm models in Mocca.
Isotherm models must implement:
- `compute_equilibrium(model, concentration, temperature)`: returns equilibrium loading
- `compute_ki(model, concentration, qstar)`: returns mass transfer coefficient
"""
abstract type AbstractIsothermModel end

"""
    compute_equilibrium(model::AbstractIsothermModel, concentration, temperature)

Compute equilibrium loading for the given concentrations and temperature.

# Arguments
- `model`: The isotherm model
- `concentration`: Vector of component concentrations [mol/m³]
- `temperature`: Temperature [K]

# Returns
- Equilibrium loading vector [mol/kg] as SVector
"""
function compute_equilibrium end

"""
    compute_ki(model::AbstractIsothermModel, concentration, qstar)

Compute mass transfer coefficients for the Linear Driving Force (LDF) model.

# Arguments
- `model`: The isotherm model
- `concentration`: Vector of component concentrations [mol/m³]
- `qstar`: Equilibrium loading vector [mol/kg]

# Returns
- Mass transfer coefficient vector [1/s] as SVector
"""
function compute_ki end

# Include specific implementations
include("haghpanah_dual_site.jl")
include("langmuir_backend.jl")
