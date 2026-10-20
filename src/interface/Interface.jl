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
using ..Framework, ..Physics, ..Balances, ..Devices, ..FixedBeds, ..Networks, ..Operation
using ..Framework: _is_solver_failure
using ..Operation: _no_targets, _plural, _simulate_schedule, _stage_forces, _sci, cycle_progress, finish_phase!, finish_progress!, plan_phases!,
    report_cycle!, report_note!, start_phase!

export AbstractIsothermData, AbstractStep, add_bed!, add_device!, add_variable!, Adsorb,
    Adsorbent, adsorbent_to_dict, attachments, axial_diffusivity, bed_dispersion, beds,
    BedSpec, Blowdown, boundary, build_isotherm, Column, compile_cycle, CompiledCycle,
    Component, component, composition, Concurrent, Cycle, cycle_schedule, cycle_steps,
    cycle_time, cyclic_steady_state, default_conductance, default_step_name, DesignVariable,
    device_role, DEVICE_ROLES, devices_at, devices_with_role, dual_adsorbent_vsa_process,
    dual_site_langmuir, DualSiteLangmuirData, end_pressure, end_pressures, Equalise,
    equilibrium_loading, Evacuate, Feed, FeedPlace, find_device, fixed_bed_model, Flowsheet,
    four_stage_vsa_process, Gas, gas_density, gas_viscosity, haghpanah_2013_column,
    haghpanah_2013_gas, Idle, initial_state, interstitial_flow_rate, isotherm,
    isotherm_table, krishnamurthy_2014_column, krishnamurthy_2014_gas, langmuir,
    library_adsorbents, library_processes, light_product_composition, load_adsorbent,
    lpp_vsa_process, mass_transfer, MoccaInputError, optimise!, ParameterPlace, performance,
    performance_gradient, Pressurise, print_flowsheet_diagram, Process, process_to_dict,
    ProcessOptions, ProcessResult, Purge, ramp, read_process, register_component!,
    register_isotherm!, register_step!, reset_step_pattern!, save_adsorbent, set_variables!,
    stages, start_pressure, step_beds, step_settings, StepContext, StepPlace,
    superficial_flow_rate, TimeStepping, two_bed_vsa_process, validate, valve,
    variable_names, variable_values, Wall, WarmStart, WET_FLUE_GAS_FEED, write_process,
    write_results

include("errors.jl")
include("gas.jl")
include("materials.jl")
include("builder.jl")
include("steps.jl")
include("process.jl")
include("variables.jl")
include("files.jl")
include("results.jl")
include("validation.jl")
include("presets.jl")

end
