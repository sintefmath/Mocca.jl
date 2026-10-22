"""
    Mocca.Physics

Material models: isotherms, mass transfer and gas thermodynamics.
"""
module Physics

import Jutul
using StaticArrays
using ..Framework

export AbstractGasThermodynamics, AbstractIsotherm, AbstractMassTransfer,
    component_enthalpies, component_heat_capacities, compute_enthalpy, compute_equilibrium,
    compute_mass_transfer_rate, DualSiteLangmuir, heat_capacity_ratio, IdealGasConstantCp,
    IdealGasShomate, LinearDrivingForce, molar_concentration, molar_enthalpy,
    molar_heat_capacity, saturation_capacity, shomate_coefficients,
    THERMODYNAMIC_REFERENCE_TEMPERATURE

include("isotherms/isotherms.jl")
include("isotherms/dual_site_langmuir.jl")
include("mass_transfer/mass_transfer.jl")
include("mass_transfer/ldf.jl")
include("thermodynamics/thermodynamics.jl")
include("thermodynamics/ideal_gas.jl")

end
