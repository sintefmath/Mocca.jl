__precompile__(true)

module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export MoccaSystem, MoccaModel, DistributedUnit, SorbentBed, FlowChannel, LumpedUnit, FlowElement
export Unit
export FixedBed, FixedBedModel, AdsorptionSystem, AdsorptionModel
export AbstractThermalModel, WithWall, Adiabatic, Isothermal
export AbstractSourceTerm
export FlowDevice, setup_device_model, AbstractDeviceLaw, Closed, LinearValve, VolumetricFlow
export PortCondition, ExponentialRamp, port_conductance
export Flowsheet, add_unit!, connect!, set_boundary!, setup_flowsheet_model
export setup_flowsheet_state, setup_flowsheet_parameters
export Stage, setup_schedule, four_stage_vsa_flowsheet, two_bed_vsa_flowsheet
export cycle_window, stream_totals, purity, recovery, productivity, sorbent_mass
export vacuum_pump_energy, specific_energy_kwh_per_tonne, component_inventory
export restart_state, cycle_change, simulate_to_cyclic_steady_state, simulate_until_mass_balance
export StreamObjective, VacuumPumpObjective, CompressorObjective, force_gradients, newton_cyclic_steady_state, cyclic_steady_state_gradient
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

export plot_state, plot_cell, plot_cell_comparison, plot_state_comparison, plot_flowsheet_streams, plot_profiles

export AbstractIsotherm, compute_equilibrium, compute_enthalpy, DualSiteLangmuir

export AbstractMassTransfer, compute_mass_transfer_rate, LinearDrivingForce
export compute_permeability, compute_dispersion

# User interface
export Units, MoccaInputError, validate
export Component, component, register_component!, Gas, Feed, composition
export AbstractIsothermData, dual_site_langmuir, langmuir, register_isotherm!
export Adsorbent, Column, Wall, equilibrium_loading, isotherm_table
export add_bed!, add_device!, devices_with_role, print_flowsheet_diagram
export AbstractStep, StepContext, step_settings, end_pressure, register_step!
export Pressurise, Adsorb, Blowdown, Evacuate, Purge, Equalise, Idle, Concurrent, Cycle, TimeStepping
export Process, simulate, cyclic_steady_state, performance, initial_state
export StepPlace, FeedPlace, ParameterPlace, add_variable!, performance_gradient, optimise!, reset_step_pattern!
export read_process, write_process, load_adsorbent, save_adsorbent, library_adsorbents, library_processes
export four_stage_vsa_process, two_bed_vsa_process, lpp_vsa_process, dual_adsorbent_vsa_process, convert_constants

import Jutul
import Jutul: simulate
import JSON
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
include("core/variables.jl")
include("core/conservation.jl")
include("core/thermal.jl")

# Input structs
include("io/constants.jl")

# L0: physics functions dispatched on physics objects
include("physics/isotherms/isotherms.jl")
include("physics/mass_transfer/mass_transfer.jl")

# L2: unit operations
include("units/fixed_bed/system.jl")
include("units/fixed_bed/setup.jl")

# L1: reusable blocks of variables, fluxes and source terms
include("blocks/geometry.jl")
include("blocks/unit_parameters.jl")
include("blocks/gas.jl")
include("blocks/sorbent.jl")
include("blocks/energy.jl")
include("blocks/wall.jl")

include("units/fixed_bed/legacy_bcs/forces.jl")
include("units/fixed_bed/select.jl")
include("units/equipment/device.jl")

# L3/L4: coupling units into flowsheets
include("coupling/ports.jl")
include("flowsheet/flowsheet.jl")

# L5: process operation
include("process/schedule.jl")
include("process/vsa.jl")
include("process/wet_flue_gas.jl")
include("process/metrics.jl")
include("process/css.jl")
include("process/gradients.jl")

# User interface: components, materials, flowsheets with roles, step
# templates, processes, design variables and files
include("interface/errors.jl")
include("interface/units.jl")
include("interface/gas.jl")
include("interface/materials.jl")
include("interface/builder.jl")
include("interface/steps.jl")
include("interface/process.jl")
include("interface/variables.jl")
include("interface/files.jl")
include("interface/validation.jl")
include("interface/presets.jl")

include("core/convergence.jl")
include("utils.jl")
include("plot.jl")
include("io/input_output.jl")
include("../models/models.jl")
end
