__precompile__(true)

module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export AdsorptionSystem, AdsorptionModel, TwoComponentAdsorptionSystem, Column
export MoccaCase

export setup_process_model
export setup_process_simulator
export setup_process_parameters
export setup_process_state
export setup_boundary_conditions
export setup_dcb_forces
export simulate_process
export mocca_domain
export column_mesh

export plot_state, plot_cell

export AbstractIsotherm, compute_equilibrium, compute_enthalpy, DualSiteLangmuir

export AbstractMassTransfer, compute_mass_transfer_rate, LinearDrivingForce
export compute_permeability, compute_dispersion

import Jutul
import JutulDarcy
using StaticArrays

import Jutul: JutulCase

const MoccaCase = JutulCase # Convenience alias for simulation cases

const GAS_CONSTANT = 8.3144598 # J/(mol·K)

const moccaResultsDir = joinpath(@__DIR__, "..", "results")

if !isdir(moccaResultsDir)
    mkpath(moccaResultsDir)
end


include("core_types/core_types.jl")

include("init/init.jl")
include("isotherms/isotherms.jl")
include("mass_transfer/mass_transfer.jl")
include("systems/systems.jl")
include("variables/variables.jl")
include("equations/equations.jl")
include("forces/forces.jl")

include("select_variable.jl")
include("convergence.jl")
include("utils.jl")
include("plot.jl")
include("input_output/input_output.jl")
include("../models/models.jl")
end
