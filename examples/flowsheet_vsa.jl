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
fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 100)
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
        Pressure = constants.p_low,
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
# Each state holds one entry per unit. The devices report the molar flow
# through them, which gives the product streams directly.
bed_states = [s[:Bed] for s in states]
f_outlet = Mocca.plot_cell(bed_states, fs[:Bed], timesteps_out, Jutul.number_of_cells(fs[:Bed].domain))

t = cumsum(timesteps_out)
F_vacuum = [s[:V_vacuum][:MolarFlow][1] for s in states]
y_CO2_extract = [s[:V_vacuum][:InletComposition][1, 1] for s in states]
n_CO2 = sum(F_vacuum .* y_CO2_extract .* timesteps_out)
println("CO₂ recovered through the vacuum line over three cycles: $(round(n_CO2, digits = 2)) mol")
