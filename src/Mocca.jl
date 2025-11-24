__precompile__(true)

module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export AdsorptionSystem, AdsorptionModel, TwoComponentAdsorptionSystem
export FlexibleAdsorptionSystem, FlexibleAdsorptionSystemWithLangmuir, FlexibleAdsorptionSystemWithHaghpanah
export AbstractIsothermModel, HaghpanahDualSiteIsotherm, LangmuirIsothermModel
export create_haghpanah_equivalent_langmuir
export MoccaCase

export setup_process_simulator
export setup_process_model
export setup_process_parameters
export setup_process_state
export setup_forces
export setup_dcb_forces
export simulate_process

export plot_state, plot_cell

import Jutul
import JutulDarcy
import CommonSolve
using StaticArrays

import Jutul: JutulCase

const MoccaCase = JutulCase # Convenience alias for simulation cases

# TODO: Remove these when n-component systems are implemented
const CO2INDEX = 1 # TODO: We don't really need this
const N2INDEX = 2 # TODO: We don't really need this

const moccaResultsDir = joinpath(@__DIR__, "..", "results")

if !isdir(moccaResultsDir)
    mkpath(moccaResultsDir)
end

# Runtime initialization function
function __init__()
    # Workaround for Langmuir.jl v0.1.2 bug where step! is not properly imported
    # This needs to happen at runtime to avoid breaking precompilation
    try
        @eval Langmuir begin
            import CommonSolve: step!
        end
    catch e
        # If this fails, it might already be imported or the issue is fixed in newer versions
        @warn "Could not inject step! into Langmuir module" exception=e
    end
end


include("core_types/core_types.jl")

include("init/init.jl")

include("isotherms/isotherms.jl")

include("systems/systems.jl")
include("variables/variables.jl")
include("equations/equations.jl")
include("forces/forces.jl")

include("select_variable.jl")
include("updates.jl")
include("updates_flexible.jl")
include("convergence.jl")
include("utils.jl")
include("plot.jl")
include("input_output/input_output.jl")
include("../models/models.jl")
end
