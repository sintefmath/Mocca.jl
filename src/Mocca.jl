__precompile__(true)

module Mocca

export initialise_state_AdsorptionColumn
export ConstantsStruct, HaghpanahConstants

export compute_permeability, calc_dispersion


export AdsorptionSystem, AdsorptionModel, TwoComponentAdsorptionSystem
export compute_equilibrium, compute_ki, initialize_from_matlab, plot_states



import Jutul
import JutulDarcy
using StaticArrays

# TODO: Remove these when n-component systems are implemented
const CO2INDEX = 1 # TODO: We don't really need this
const N2INDEX = 2 # TODO: We don't really need this



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
