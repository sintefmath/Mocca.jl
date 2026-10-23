module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export MoccaSystem, MoccaModel, DistributedUnit, SorbentBed, FlowChannel, LumpedUnit, FlowElement
export Unit
export FixedBed, FixedBedModel, AdsorptionSystem, AdsorptionModel
export AbstractThermalModel, WithWall, Adiabatic, Isothermal
export AbstractSourceTerm
export FlowDevice, setup_device_model, AbstractDeviceLaw, Closed, LinearValve, VolumetricFlow
export PortCondition, ExponentialRamp, port_conductance
export Flowsheet, UnitNetwork, add_unit!, connect!, set_boundary!, setup_flowsheet_model
export AbstractConnection, StreamConnection, PortRef
export setup_flowsheet_state, setup_flowsheet_parameters
export Stage, setup_schedule
export cycle_window, stream_totals, purity, recovery, productivity, adsorbent_volume, sorbent_mass
export vacuum_pump_energy, specific_energy_kwh_per_tonne, component_inventory
export restart_state, cycle_change, simulate_to_cyclic_steady_state, simulate_until_mass_balance, CycleProgress
export StreamObjective, VacuumPumpObjective, CompressorObjective, force_gradients, newton_cyclic_steady_state, cyclic_steady_state_gradient
export number_of_components, component_names
export MoccaCase

export setup_process_model
export setup_process_simulator
export setup_process_parameters
export setup_process_state
export simulate_process
export mocca_domain
export column_mesh

export plot_state, plot_cell, plot_cell_comparison, plot_state_comparison, plot_flowsheet_streams, plot_profiles

export AbstractIsotherm, compute_equilibrium, compute_enthalpy, DualSiteLangmuir

export AbstractMassTransfer, compute_mass_transfer_rate, LinearDrivingForce

export AbstractGasThermodynamics, IdealGasConstantCp, IdealGasShomate, component_heat_capacities,
    component_enthalpies, molar_heat_capacity, molar_enthalpy, molar_concentration, heat_capacity_ratio
export AbstractGasTransport, ConstantTransport, SutherlandTransport, sutherland_transport, TransportReference, mixture_viscosity,
    mixture_thermal_conductivity, diffusivity_scale
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
export read_process, write_process, write_results, load_adsorbent, save_adsorbent, library_adsorbents, library_processes
export four_stage_vsa_process, two_bed_vsa_process, lpp_vsa_process, dual_adsorbent_vsa_process, convert_constants

export haghpanah_DCB_input, haghpanah_cyclic_input

import Jutul
import Jutul: simulate

# Each module only uses the modules above it. A module that adds a method to a
# function of another writes the function's full name, Framework.ports(...).
include("Units.jl")
include("framework/Framework.jl")
include("physics/Physics.jl")
include("balances/Balances.jl")
include("devices/Devices.jl")
include("fixedbeds/FixedBeds.jl")
include("networks/Networks.jl")
include("operation/Operation.jl")
include("interface/Interface.jl")
include("legacy/Legacy.jl")
include("plotting/Plotting.jl")

using .Framework, .Physics, .Balances, .Devices, .FixedBeds, .Networks, .Operation
using .Interface, .Legacy, .Plotting

end
