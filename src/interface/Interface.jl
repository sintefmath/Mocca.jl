"""
    Mocca.Interface

The user interface: components, materials, flowsheets with roles, steps,
processes, design variables and files.
"""
module Interface

import Jutul
import Jutul: simulate
using StaticArrays
import JSON
import ..Units
using ..Framework, ..Physics, ..Balances, ..Devices, ..FixedBeds, ..MovingBeds, ..Networks, ..Operation
using ..Framework: _is_solver_failure
using ..Operation: _no_targets, _plural, _simulate_schedule, _stage_forces, _sci, cycle_progress, finish_phase!, finish_progress!, plan_phases!,
    report_cycle!, report_note!, start_phase!, DEFAULT_NONLINEAR_TOLERANCE, STEADY_STATE_NONLINEAR_TOLERANCE

export AbstractIsothermData, AbstractStep, add_bed!, add_device!, add_integral!,
    add_measure!, add_moving_bed!, add_variable!, Adsorb, Adsorbent, adsorbent_to_dict,
    attachments, axial_diffusivity, bed_dispersion, beds, BedSpec, Blowdown, boundary,
    build_isotherm, chapman_enskog_transport, check_gradients, Column, compile_cycle,
    CompiledCycle, component, Component, composition, Concurrent, Cycle, cycle_schedule,
    cycle_steps, cycle_time, cyclic_steady_state, CyclicProcess, CyclicResult,
    default_conductance, default_step_name, DesignVariable, device_role, DEVICE_ROLES,
    devices_at, devices_with_role, dual_adsorbent_vsa_process, dual_site_langmuir,
    DualSiteLangmuirData, end_pressure, end_pressures, Equalise, equilibrium_loading,
    Evacuate, explore, Feed, FeedPlace, find_device, fixed_bed_model, Flowsheet,
    four_stage_vsa_process, Gas, gas_density, gas_diffusivities, gas_heat_capacity_ratio,
    gas_viscosity, GeometryPlace, haghpanah_2013_column, haghpanah_2013_gas, Idle,
    initial_state, interstitial_flow_rate, isotherm, isotherm_table,
    krishnamurthy_2014_column, krishnamurthy_2014_gas, langmuir, library_adsorbents,
    library_processes, light_product_composition, load_adsorbent, lpp_vsa_process,
    mass_transfer, MoccaInputError, moving_bed_model, MovingBedSpec, operating_stage,
    optimise!, ParameterPlace, pareto_front!, performance, performance_gradient,
    performance_integrals, performance_measures, performance_quantities, performance_sums,
    Place, Pressurise, print_flowsheet_diagram, Process, process_to_dict,
    pump_heat_capacity_ratio,
    ProcessOptions, ProcessResult, Purge, ramp, read_process, register_component!,
    register_isotherm!, register_step!, reset_step_pattern!, save_adsorbent, ScaledPlace, set_variables!, SettingPlace,
    stages, start_pressure, steady_state, SteadyProcess, SteadyResult, step_beds, step_settings,
    StepContext, StepPlace, superficial_flow_rate, TimeStepping, two_bed_vsa_process, validate,
    valve, variable_names, variable_values, virial, VirialData, Wall, WarmStart, WET_FLUE_GAS_FEED, write_process,
    write_results

include("errors.jl")
include("gas.jl")
include("materials.jl")
include("builder.jl")
include("steps.jl")
include("process.jl")
include("variables.jl")
include("steady.jl")
include("measures.jl")
include("optimisation.jl")
include("checks.jl")
include("files.jl")
include("results.jl")
include("validation.jl")
include("presets.jl")

end
