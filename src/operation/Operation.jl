"""
    Mocca.Operation

Running a flowsheet: stages, cyclic steady state, performance metrics and
adjoint gradients.
"""
module Operation

import Jutul
using LinearAlgebra
using ..Framework, ..Physics, ..Devices, ..Networks

export adsorbent_volume, component_inventory, compressor_energy, CompressorObjective,
    cycle_change, CycleProgress, cycle_window, cyclic_steady_state_gradient, force_gradients,
    LOADING_SCALE_FRACTION, MOLE_FRACTION_SCALE_FLOOR, newton_cyclic_steady_state,
    parallel_stages, productivity, purity, recovery, restart_state, setup_process_simulator,
    setup_schedule, simulate_process, simulate_to_cyclic_steady_state,
    simulate_until_mass_balance, sorbent_mass, specific_energy_kwh_per_tonne, Stage,
    stream_totals, StreamObjective, vacuum_pump_energy, VacuumPumpObjective

include("simulator.jl")
include("schedule.jl")
include("metrics.jl")
include("progress.jl")
include("css.jl")
include("gradients.jl")

end
