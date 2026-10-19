#=
# Cyclic VSA as a flowsheet

This example looks inside the four-stage vacuum swing adsorption (VSA) cycle
of Haghpanah et al. (2013) as a flowsheet of connected units, and compares it
with the stage-specific boundary conditions of earlier Mocca versions. A fixed
bed is connected to three flow devices:

    feed ──V_feed──┐                 ┌──V_product── raffinate
                   bottom  Bed  top
    extract ──V_vacuum──┘

Each stage of the cycle is a set of device settings (valve open or closed,
flow controller rate, boundary pressure profile). Nothing inside the bed
changes between stages, so the same bed model can be wired into larger
flowsheets, such as several beds sharing equalisation lines.
=#

import Mocca
using Mocca.Units

# # 1. The process
#
# `four_stage_vsa_process` builds the flowsheet and its cycle of steps; the
# [cyclic VSA](cyclic_vsa_haghpanah_2013_co2_n2.md) example builds the same
# from its parts. The grid, stage durations, 1 s time steps and initial state
# (nitrogen at 1 bar) are those of the boundary-condition version, so the two
# can be compared at the end.
ncells = 200
process = Mocca.four_stage_vsa_process(ncells = ncells, timestepping = Mocca.TimeStepping(max_dt = 1.0, first_dt = 1.0),
    initial_state = Dict(:Bed => (pressure = 1bar,)))

# The steps compile into stages, each a set of device settings:
for st in Mocca.stages(process)
    println(st.name, " (", st.duration, " s): ")
    for (d, s) in st.settings
        println("    ", d, " => ", s)
    end
end

# # 2. Simulate
result = Mocca.simulate(process; cycles = 3, info_level = 0)

# # 3. Look at the results
#
# Each state holds one entry per unit. The bed states can be plotted with the
# same functions as a single-column simulation. First the primary variables at
# the outlet (top) of the bed through time:
fs = process.flowsheet
bed = fs[:Bed]
bed_states = [s[:Bed] for s in result.states]
f_outlet = Mocca.plot_cell(bed_states, bed, result.timesteps, ncells)

# and along the bed at the end of the simulation:
f_column = Mocca.plot_state(bed_states[end], bed)

# The devices report the molar flow through them and the composition at their
# ports, which gives the streams entering and leaving the bed. Flow is positive
# from a device's inlet to its outlet, so the vacuum line is positive when gas
# is drawn out of the bed.
f_streams = Mocca.plot_flowsheet_streams(result.states, process.model, result.timesteps; stages = Mocca.stages(process))

# The stream totals give the CO₂ collected through the vacuum line.
n_extract = Mocca.stream_totals(result.states, result.timesteps, :V_vacuum)
println("CO₂ recovered through the vacuum line over three cycles: $(round(n_extract[1], digits = 2)) mol")

# # 4. Comparison with the boundary-condition version
#
# Run the same cycle with the stage-specific boundary conditions of earlier
# Mocca versions. They are deprecated, and kept for comparisons like this.
constants = Mocca.HaghpanahConstants{Float64}()
bc_model = Mocca.setup_process_model(Mocca.FixedBed(constants), constants; ncells = ncells)
bc_state0 = Mocca.setup_process_state(bc_model;
    Pressure = 1bar,
    Temperature = 298.15,
    WallTemperature = constants.T_a,
    y = [1e-10, 1.0 - 1e-10],
)
bcs = Mocca.setup_boundary_conditions(constants, ["pressurisation", "adsorption", "blowdown", "evacuation"])
bc_forces, bc_timesteps = Mocca.setup_forces(bc_model, [st.duration for st in Mocca.stages(process)], bcs; num_cycles = 3, max_dt = 1)
bc_case = Mocca.MoccaCase(bc_model, bc_timesteps, bc_forces; state0 = bc_state0, parameters = Mocca.setup_process_parameters(bc_model))
bc_states, bc_timesteps_out = Mocca.simulate_process(bc_case; info_level = 0, output_substates = true)

runs = [
    "Flowsheet" => (bed_states, result.timesteps),
    "Boundary conditions" => (bc_states, bc_timesteps_out),
]

# At the outlet and at the feed end the two agree closely. CO₂ has not reached
# the outlet within three cycles. The small differences that remain come from
# the time stepping: the solver splits the first step of some stages
# differently in the two simulations.
f_outlet_comparison = Mocca.plot_cell_comparison(runs, bed, ncells)
f_inlet_comparison = Mocca.plot_cell_comparison(runs, bed, 1)

# Profiles along the bed at the end of the simulation:
f_column_comparison = Mocca.plot_state_comparison(
    [label => first(r)[end] for (label, r) in runs], bed)
