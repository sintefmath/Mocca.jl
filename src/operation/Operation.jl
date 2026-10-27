"""
    Mocca.Operation

Running a flowsheet: stages, cyclic steady state, performance metrics and
adjoint gradients.
"""
module Operation

import Jutul
using LinearAlgebra
using ..Framework, ..Physics, ..Devices, ..Networks

export adsorbent_volume, blower_energy, BlowerObjective, component_inventory, compressor_energy, CompressorObjective,
    cycle_change, CycleProgress, cycle_window, cyclic_steady_state_gradient, equivalent_work, EquivalentWorkObjective,
    force_gradients, heat_duty, HeatDutyObjective, hx_fluid_inlet_temperature, integral_values, LOADING_SCALE_FRACTION,
    measure_partials, MOLE_FRACTION_SCALE_FLOOR, moving_sorbent_units, newton_cyclic_steady_state, objective_state,
    parallel_stages, productivity, purity, recovery, restart_state, setup_process_simulator, setup_schedule,
    simulate_process, simulate_to_cyclic_steady_state, simulate_until_mass_balance, sorbent_circulation_time,
    sorbent_inventory, sorbent_mass, specific_energy_kwh_per_tonne, Stage, stream_rates, stream_totals,
    StreamObjective, SumObjective, vacuum_pump_energy, VacuumPumpObjective

include("simulator.jl")
include("schedule.jl")
include("metrics.jl")
include("progress.jl")
include("css.jl")
include("gradients.jl")
include("measures.jl")

end
