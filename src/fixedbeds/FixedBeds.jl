"""
    Mocca.FixedBeds

The fixed bed: a column of stationary sorbent with gas flowing through it.
"""
module FixedBeds

import Jutul
using StaticArrays
using ..Framework, ..Physics, ..Balances
using ..Balances: _is_axisymmetric

export AdsorptionModel, AdsorptionSystem, compute_dispersion, compute_permeability,
    compute_radial_dispersion, FixedBed, FixedBedModel

include("system.jl")
include("correlations.jl")
include("select.jl")

end
