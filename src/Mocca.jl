__precompile__(true)

module Mocca
export AdsorptionFlowSystem, AdsorptionFlowModel, compute_equilibrium, compute_ki, initialize_from_matlab, plot_states
import Jutul
import JutulDarcy
using StaticArrays, ForwardDiff

const CO2INDEX = 1 # TODO: We don't really need this
const N2INDEX = 2 # TODO: We don't really need this
include("parameters.jl")
include("system.jl")
include("variable_structs.jl")
include("select_variable.jl")
include("updates.jl")
include("flux/flux.jl")
include("convergence.jl")
include("bc.jl")
include("bc_adsorption.jl")
include("init.jl")
include("plot.jl")
end
