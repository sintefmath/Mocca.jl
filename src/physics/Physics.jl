"""
    Mocca.Physics

Material models: isotherms, mass transfer, gas thermodynamics and gas
transport properties.
"""
module Physics

import Jutul
using StaticArrays
using ..Framework

export AbstractGasThermodynamics, AbstractGasTransport, AbstractIsotherm,
    AbstractMassTransfer, component_enthalpies, component_heat_capacities,
    component_thermal_conductivities, component_viscosities, compute_enthalpy,
    compute_equilibrium, compute_mass_transfer_rate, conductivity_factor,
    ConstantTransport, diffusivity_factor, diffusivity_scale, DualSiteLangmuir,
    heat_capacity_ratio, IdealGasConstantCp, IdealGasShomate, LinearDrivingForce,
    mixture_thermal_conductivity, mixture_viscosity, molar_concentration, molar_enthalpy,
    molar_heat_capacity, saturation_capacity, shomate_coefficients, SutherlandTransport,
    is_constant_transport, THERMODYNAMIC_REFERENCE_TEMPERATURE, TransportReference,
    viscosity_factor

include("isotherms/isotherms.jl")
include("isotherms/dual_site_langmuir.jl")
include("mass_transfer/mass_transfer.jl")
include("mass_transfer/ldf.jl")
include("thermodynamics/thermodynamics.jl")
include("thermodynamics/ideal_gas.jl")
include("transport/transport.jl")
include("transport/sutherland.jl")

end
