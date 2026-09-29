__precompile__(true)

module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export MoccaSystem, MoccaModel, DistributedUnit, SorbentBed, FlowChannel, LumpedUnit, FlowElement
export Unit, Column
export FixedBed, FixedBedModel, AdsorptionSystem, AdsorptionModel
export number_of_components, component_names
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
using StaticArrays

import Jutul: JutulCase

const MoccaCase = JutulCase # Convenience alias for simulation cases

const GAS_CONSTANT = 8.3144598 # J/(mol·K)

const moccaResultsDir = joinpath(@__DIR__, "..", "results")

if !isdir(moccaResultsDir)
    mkpath(moccaResultsDir)
end


# Core: system hierarchy, entities, shared interfaces
include("core/types.jl")

# Input structs
include("io/constants.jl")

# L0: physics functions dispatched on physics objects
include("physics/isotherms/isotherms.jl")
include("physics/mass_transfer/mass_transfer.jl")

# L2: unit operations
include("units/fixed_bed/system.jl")
include("units/fixed_bed/setup.jl")

# Variables and equations (to be split into reusable blocks)
include("variables/variables.jl")
include("equations/equations.jl")
include("units/fixed_bed/legacy_bcs/forces.jl")
include("units/fixed_bed/select.jl")

include("core/convergence.jl")
include("utils.jl")
include("plot.jl")
include("io/input_output.jl")
include("../models/models.jl")
end
