"""
    Mocca.FixedBeds

The fixed bed: a column of stationary sorbent with gas flowing through it.
"""
module FixedBeds

import Jutul
using StaticArrays
using ..Framework, ..Physics, ..Balances
using ..Balances: _thermodynamics

export AdsorptionModel, AdsorptionSystem, FixedBed, FixedBedModel

include("system.jl")
include("select.jl")

end
