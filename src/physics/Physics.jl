"""
    Mocca.Physics

Material models: isotherms, mass transfer, gas thermodynamics, gas transport
properties and axial dispersion.
"""
module Physics

import Jutul
using StaticArrays
using ..Framework

export AbstractAxialDispersion, AbstractGasThermodynamics, AbstractGasTransport,
    AbstractIsotherm, AbstractMassTransfer, binary_diffusivity, ChapmanEnskogTransport,
    component_enthalpies, component_heat_capacities, component_viscosities, compute_enthalpy,
    compute_equilibrium, compute_mass_transfer_rate, ConstantDispersion, ConstantTransport,
    diffusivity_factors, dispersion_coefficient, DualSiteLangmuir, gab_loading, has_knudsen,
    heat_capacity_ratio, IdealGasConstantCp, IdealGasPolynomial, IdealGasShomate,
    is_constant_transport, LinearDrivingForce, LumpedLinearDrivingForce,
    minimum_mole_fraction, mixture_viscosity, molar_concentration, molar_enthalpy,
    molar_heat_capacity, molecular_diffusivities, polynomial_coefficients, pore_diffusivity,
    saturation_capacity, shomate_coefficients, THERMODYNAMIC_REFERENCE_TEMPERATURE, TothGAB,
    TransportReference, Virial, viscosity_factor, WakaoFunazkri, water_vapour_pressure

include("isotherms/isotherms.jl")
include("isotherms/dual_site_langmuir.jl")
include("isotherms/virial.jl")
include("isotherms/toth_gab.jl")
include("mass_transfer/mass_transfer.jl")
include("mass_transfer/ldf.jl")
include("thermodynamics/thermodynamics.jl")
include("thermodynamics/ideal_gas.jl")
include("transport/transport.jl")
include("transport/chapman_enskog.jl")
include("transport/dispersion.jl")

end
