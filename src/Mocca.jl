__precompile__(true)

module Mocca

export ConstantsStruct, HaghpanahConstants
export AdsorptionSystem, AdsorptionModel, TwoComponentAdsorptionSystem
export MoccaCase

export setup_adsorption_simulator
export setup_adsorption_model
export setup_adsorption_parameters
export setup_adsorption_state
export simulate_adsorption

export plot_state, plot_cell

import Jutul
import JutulDarcy
using StaticArrays

import Jutul: JutulCase

const MoccaCase = JutulCase # Convenience alias for simulation cases

# TODO: Remove these when n-component systems are implemented
const CO2INDEX = 1 # TODO: We don't really need this
const N2INDEX = 2 # TODO: We don't really need this


include("core_types/core_types.jl")

include("init/init.jl")
include("systems/systems.jl")
include("variables/variables.jl")
include("equations/equations.jl")
include("forces/forces.jl")

include("select_variable.jl")
include("updates.jl")
include("convergence.jl")
include("utils.jl")
include("plot.jl")
end
