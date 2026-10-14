"""
    Mocca.FixedBeds

The fixed bed: a column of stationary sorbent with gas flowing through it.
"""
module FixedBeds

import Jutul
using StaticArrays
using ..Framework, ..Physics, ..Balances

export AdsorptionModel, AdsorptionSystem, compute_dispersion, compute_permeability, FixedBed,
    FixedBedModel

include("system.jl")
include("correlations.jl")
include("select.jl")

end
