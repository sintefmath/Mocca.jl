module Mocca

export ConstantsStruct, HaghpanahConstants, InfoStruct
export MoccaSystem, MoccaModel, DistributedUnit, SorbentBed, FlowChannel, LumpedUnit, FlowElement
export Unit
export FixedBed, FixedBedModel, AdsorptionSystem, AdsorptionModel
export MovingBed, MovingBedModel, setup_moving_bed_model, tube_bundle, moving_bed_sorbent_htc
export AbstractLimiter, limited_slope, Upwind, VanAlbada, VanLeer, Minmod
export AbstractThermalModel, WithWall, Adiabatic, Isothermal, HeatExchanger
export AbstractSourceTerm
export FlowDevice, setup_device_model, AbstractDeviceLaw, Closed, LinearValve, VolumetricFlow
export FixedMolarFlow, BackPressureController
export PortCondition, ExponentialRamp, port_conductance
export Flowsheet, UnitNetwork, add_unit!, connect!, set_boundary!, setup_flowsheet_model
export AbstractConnection, StreamConnection, PortRef
export SorbentConnection, HeatTransferFluidConnection, connect_solids!, connect_heat_transfer_fluid!, sorbent_path
export setup_flowsheet_state, setup_flowsheet_parameters
export Stage, setup_schedule
export cycle_window, stream_totals, purity, recovery, productivity, adsorbent_volume, sorbent_mass
export stream_rates, heat_duty, equivalent_work, sorbent_inventory, sorbent_circulation_time
export vacuum_pump_energy, compressor_energy, blower_energy, specific_energy_kwh_per_tonne, component_inventory
export restart_state, cycle_change, simulate_to_cyclic_steady_state, simulate_until_mass_balance, CycleProgress
export StreamObjective, VacuumPumpObjective, CompressorObjective, BlowerObjective, HeatDutyObjective, JacketHeatObjective, EquivalentWorkObjective
export SumObjective, objective_state, integral_values, measure_partials
export force_gradients, newton_cyclic_steady_state, cyclic_steady_state_gradient
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
export plot_moving_bed_profiles

export AbstractIsotherm, compute_equilibrium, compute_enthalpy, DualSiteLangmuir, Virial, TothGAB, water_vapour_pressure, minimum_mole_fraction

export AbstractMassTransfer, compute_mass_transfer_rate, LinearDrivingForce, LumpedLinearDrivingForce

export AbstractGasThermodynamics, IdealGasConstantCp, IdealGasShomate, component_heat_capacities,
    component_enthalpies, molar_heat_capacity, molar_enthalpy, molar_concentration, heat_capacity_ratio
export AbstractGasTransport, ConstantTransport, ChapmanEnskogTransport, chapman_enskog_transport, TransportReference,
    mixture_viscosity, molecular_diffusivities, IdealGasPolynomial, AbstractAxialDispersion, ConstantDispersion, WakaoFunazkri
export compute_permeability, compute_dispersion

# User interface
export Units, MoccaInputError, validate
export Component, component, register_component!, Gas, Feed, composition
export AbstractIsothermData, dual_site_langmuir, langmuir, virial, toth_gab, lumped_ldf, register_isotherm!
export Adsorbent, Column, Wall, equilibrium_loading, isotherm_table
export add_bed!, add_moving_bed!, add_device!, devices_with_role, print_flowsheet_diagram
export AbstractStep, StepContext, step_settings, end_pressure, register_step!
export Pressurise, Adsorb, Blowdown, Evacuate, Purge, Equalise, Idle, Concurrent, Cycle, TimeStepping
export Process, CyclicProcess, ProcessResult, CyclicResult, simulate, cyclic_steady_state, performance, jacket_heat, initial_state
export SteadyProcess, SteadyResult, steady_state
export Place, StepPlace, FeedPlace, ParameterPlace, SettingPlace, GeometryPlace, ScaledPlace, add_variable!, performance_gradient, optimise!, pareto_front!, explore, check_gradients, reset_step_pattern!
export performance_integrals, performance_quantities, performance_measures, performance_sums, add_integral!, add_measure!
export read_process, write_process, write_results, load_adsorbent, save_adsorbent, library_adsorbents, library_processes
export four_stage_vsa_process, two_bed_vsa_process, lpp_vsa_process, dual_adsorbent_vsa_process, dac_tvsa_process, convert_constants

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
include("movingbeds/MovingBeds.jl")
include("networks/Networks.jl")
include("operation/Operation.jl")
include("interface/Interface.jl")
include("legacy/Legacy.jl")
include("plotting/Plotting.jl")

using .Framework, .Physics, .Balances, .Devices, .FixedBeds, .MovingBeds, .Networks, .Operation
using .Interface, .Legacy, .Plotting

end
