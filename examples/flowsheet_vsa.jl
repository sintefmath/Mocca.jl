#=
# Cyclic VSA as a flowsheet

This example builds the four-stage vacuum swing adsorption (VSA) cycle of
Haghpanah et al. (2013) from connected units instead of stage-specific
boundary conditions. A fixed bed is connected to three flow devices:

    feed ──V_feed──┐                 ┌──V_product── raffinate
                   bottom  Bed  top
    extract ──V_vacuum──┘

Each stage of the cycle is a set of device settings (valve open or closed,
flow controller rate, boundary pressure profile). Nothing inside the bed
changes between stages, so the same bed model can be wired into larger
flowsheets, such as several beds sharing equalisation lines.
=#

import Jutul
import Jutul: si_unit
import Mocca

constants = Mocca.HaghpanahConstants{Float64}()

# # 1. Build the flowsheet
#
# `four_stage_vsa_flowsheet` wires the bed and devices and returns the stages.
# The steps it performs are:
#
# ```julia
# fs = Mocca.Flowsheet()
# Mocca.add_unit!(fs, :Bed, bed_model)
# Mocca.add_unit!(fs, :V_feed, Mocca.setup_device_model(["CO2", "N2"]))
# Mocca.connect!(fs, (:Bed, :bottom) => (:V_feed, :outlet))
# Mocca.set_boundary!(fs, :V_feed, :inlet, Mocca.PortCondition(pressure = 1e5, temperature = 298.15, composition = [0.15, 0.85]))
# ```
#
# The grid, stage durations and initial state are those of the
# [cyclic VSA](cyclic_vsa_haghpanah_2013_co2_n2.md) example, so the two can be
# compared at the end.
ncells = 200
fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = ncells)
model = Mocca.setup_flowsheet_model(fs)

for st in stages
    println(st.name, " (", st.duration, " s): ", join(keys(st.settings), ", "))
end

# # 2. Initial state, parameters and schedule
#
# Only units with holdup need an initial state; device states follow from what
# their ports are attached to.
state0 = Mocca.setup_flowsheet_state(fs;
    Bed = Mocca.setup_process_state(fs[:Bed];
        Pressure = 1*si_unit(:bar),
        Temperature = 298.15,
        WallTemperature = constants.T_a,
        y = [1e-10, 1.0 - 1e-10],
    )
)
parameters = Mocca.setup_flowsheet_parameters(fs)
forces, timesteps = Mocca.setup_schedule(fs, model, stages; num_cycles = 3, max_dt = 1.0)

# # 3. Simulate
case = Mocca.MoccaCase(model, timesteps, forces; state0 = state0, parameters = parameters)
states, timesteps_out = Mocca.simulate_process(case; info_level = 0, output_substates = true)

# # 4. Look at the results
#
# Each state holds one entry per unit. The bed states can be plotted with the
# same functions as a single-column simulation. First the primary variables at
# the outlet (top) of the bed through time:
bed = fs[:Bed]
bed_states = [s[:Bed] for s in states]
f_outlet = Mocca.plot_cell(bed_states, bed, timesteps_out, ncells)

# and along the bed at the end of the simulation:
f_column = Mocca.plot_state(bed_states[end], bed)

# The devices report the molar flow through them and the composition at their
# ports, which gives the streams entering and leaving the bed. Flow is positive
# from a device's inlet to its outlet, so the vacuum line is positive when gas
# is drawn out of the bed.
f_streams = Mocca.plot_flowsheet_streams(states, model, timesteps_out; stages = stages)

# The stream totals give the CO₂ collected through the vacuum line.
n_extract = Mocca.stream_totals(states, timesteps_out, :V_vacuum)
println("CO₂ recovered through the vacuum line over three cycles: $(round(n_extract[1], digits = 2)) mol")

# # 5. Comparison with the boundary-condition version
#
# Run the same cycle with the stage-specific boundary conditions of the
# [cyclic VSA](cyclic_vsa_haghpanah_2013_co2_n2.md) example.
bc_model = Mocca.setup_process_model(Mocca.AdsorptionSystem(constants), constants; ncells = ncells)
bc_state0 = Mocca.setup_process_state(bc_model;
    Pressure = 1*si_unit(:bar),
    Temperature = 298.15,
    WallTemperature = constants.T_a,
    y = [1e-10, 1.0 - 1e-10],
)
bcs = Mocca.setup_boundary_conditions(constants, ["pressurisation", "adsorption", "blowdown", "evacuation"])
bc_forces, bc_timesteps = Mocca.setup_forces(bc_model, [st.duration for st in stages], bcs; num_cycles = 3, max_dt = 1)
bc_case = Mocca.MoccaCase(bc_model, bc_timesteps, bc_forces; state0 = bc_state0, parameters = Mocca.setup_process_parameters(bc_model))
bc_states, bc_timesteps_out = Mocca.simulate_process(bc_case; info_level = 0, output_substates = true)

runs = [
    "Flowsheet" => (bed_states, timesteps_out),
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
