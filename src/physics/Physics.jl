"""
    Mocca.Physics

Material models: isotherms and mass transfer.
"""
module Physics

import Jutul
using StaticArrays
using ..Framework

export AbstractIsotherm, AbstractMassTransfer, compute_enthalpy, compute_equilibrium,
    compute_mass_transfer_rate, DualSiteLangmuir, LinearDrivingForce, saturation_capacity

include("isotherms/isotherms.jl")
include("isotherms/dual_site_langmuir.jl")
include("mass_transfer/mass_transfer.jl")
include("mass_transfer/ldf.jl")

end
