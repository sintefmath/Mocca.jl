__precompile__(true)

module Mocca
export AdsorptionSystem, AdsorptionModel, TwoComponentAdsorptionSystem
export compute_equilibrium, compute_ki, initialize_from_matlab, plot_states
export compute_permeability, calc_dispersion
import Jutul
import JutulDarcy
using StaticArrays, ForwardDiff

const CO2INDEX = 1 # TODO: We don't really need this
const N2INDEX = 2 # TODO: We don't really need this
include("parameters.jl")
include("systems/AdsorptionSystem.jl")
include("systems/TwoComponentAdsorption.jl")
include("utils.jl")
include("primary_variables.jl")
include("secondary_variables.jl")
include("select_variable.jl")
include("updates.jl")
include("equations/flux.jl")
include("convergence.jl")
include("forces/bc_pressurisation.jl")
include("forces/bc_adsorption.jl")
include("forces/bc_blowdown.jl")
include("forces/bc_evacuation.jl")
include("init.jl")
include("plot.jl")
end
