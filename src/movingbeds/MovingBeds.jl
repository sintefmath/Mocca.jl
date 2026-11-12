"""
    Mocca.MovingBeds

The moving bed: a column whose sorbent moves down through it while gas flows
through, the unit of moving bed temperature swing adsorption (MBTSA).
"""
module MovingBeds

import Jutul
using StaticArrays
using ..Framework, ..Physics, ..Balances
using ..Balances: _thermodynamics
using ..Framework: _gas_state

export moving_bed_sorbent_htc, MovingBed, MovingBedModel, setup_moving_bed_model, tube_bundle

include("system.jl")
include("setup.jl")

end
