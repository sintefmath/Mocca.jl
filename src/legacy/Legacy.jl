"""
    Mocca.Legacy

Inputs of Mocca 0.1.0, and moving them over to the user interface.
"""
module Legacy

import Jutul
using StaticArrays
using JSON
using Parameters
using ..Framework, ..Interface

export adsorptionConstants, ConstantsStruct, convert_constants, export_cell_results,
    haghpanah_cyclic_input, haghpanah_DCB_input, HaghpanahConstants, InfoStruct, parse_input,
    processInfo

include("constants.jl")
include("convert_constants.jl")
include("export_cell_results.jl")
include("parse_input.jl")
include("../../models/julia/model_library.jl")

end
