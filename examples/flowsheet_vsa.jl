#=
# Cyclic VSA as a flowsheet

This example looks inside the four-stage vacuum swing adsorption (VSA) cycle
of Haghpanah et al. (2013) as a flowsheet of connected units, and then builds
the same flowsheet by hand from the layer underneath. A fixed bed is connected
to three flow devices:

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
# from its parts. The run uses steps of 1 s and starts from nitrogen at 1 bar.
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

# # 4. The same flowsheet by hand
#
# A [`Flowsheet`](@ref Mocca.Flowsheet) holds a `UnitNetwork`: the units, the
# connections between their ports and the boundary conditions of device ports
# that are not connected. Building the
# network directly means giving everything the process filled in from the
# device roles and steps: the bed model, the boundaries, the valve sizes and
# the device settings of every stage.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = [0.15, 0.85], pressure = 1bar)
bed_model = Mocca.fixed_bed_model(gas, feed, Mocca.load_adsorbent("zeolite_13x_haghpanah_2013"),
    Mocca.haghpanah_2013_column(); ncells = ncells, dispersion_velocity = 1.0)

net = Mocca.UnitNetwork()
Mocca.add_unit!(net, :Bed, bed_model)
for d in (:V_feed, :V_product, :V_vacuum)
    Mocca.add_unit!(net, d, Mocca.setup_device_model(gas.names))
end
Mocca.connect!(net, (:Bed, :bottom) => (:V_feed, :outlet))
Mocca.connect!(net, (:Bed, :top) => (:V_product, :inlet))
Mocca.connect!(net, (:Bed, :bottom) => (:V_vacuum, :inlet))

# The free ports are the feed, the raffinate line and the vacuum line:
T, y_feed, y_light = 298.15, feed.composition, feed.light_product
condition(P, y) = Mocca.PortCondition(pressure = P, temperature = T, composition = y)
Mocca.set_boundary!(net, :V_feed, :inlet, condition(1bar, y_feed))
Mocca.set_boundary!(net, :V_product, :outlet, condition(1bar, y_light))
Mocca.set_boundary!(net, :V_vacuum, :outlet, condition(0.1bar, y_feed))

# Each stage gives the settings of the devices it opens; the others stay
# closed. Valves have the conductance of the half cell at the bed end they are
# attached to, and pressures follow exponential ramps at 0.5 s⁻¹.
C_in, C_out = Mocca.port_conductance(bed_model, 1), Mocca.port_conductance(bed_model, ncells)
ramp(from, to) = Mocca.ExponentialRamp(from, to, 0.5)
rate = Mocca.interstitial_flow_rate(bed_model, 1.0)
stages = [
    Mocca.Stage("pressurisation", 15.0;
        V_feed = (law = Mocca.LinearValve(C_in), inlet = condition(ramp(0.1bar, 1bar), y_feed))),
    Mocca.Stage("adsorption", 15.0;
        V_feed = (law = Mocca.VolumetricFlow(rate), inlet = condition(1bar, y_feed)),
        V_product = (law = Mocca.LinearValve(C_out), outlet = condition(1bar, y_light))),
    Mocca.Stage("blowdown", 30.0;
        V_product = (law = Mocca.LinearValve(C_out), outlet = condition(ramp(1bar, 0.2bar), y_light))),
    Mocca.Stage("evacuation", 40.0;
        V_vacuum = (law = Mocca.LinearValve(C_in), outlet = condition(ramp(0.2bar, 0.1bar), y_feed))),
]

# The network becomes a Jutul `MultiModel`, the stages its forces, and it runs
# like any Jutul case:
model = Mocca.setup_flowsheet_model(net)
bed0 = Mocca.setup_unit_state(bed_model; pressure = 1bar, temperature = T, composition = y_light)
state0 = Mocca.setup_flowsheet_state(net; Bed = bed0)
parameters = Mocca.setup_flowsheet_parameters(net)
forces, dt = Mocca.setup_schedule(net, model, stages; num_cycles = 3, max_dt = 1.0)
states, timesteps = Mocca.simulate_process(Mocca.MoccaCase(model, dt, forces; state0 = state0, parameters = parameters);
    info_level = 0, output_substates = true)

runs = [
    "Process" => (bed_states, result.timesteps),
    "By hand" => ([s[:Bed] for s in states], timesteps),
]

# The two describe the same flowsheet and give the same results, at the
# outlet and at the feed end:
f_outlet_comparison = Mocca.plot_cell_comparison(runs, bed, ncells)
f_inlet_comparison = Mocca.plot_cell_comparison(runs, bed, 1)

# and along the bed at the end of the simulation:
f_column_comparison = Mocca.plot_state_comparison(
    [label => first(r)[end] for (label, r) in runs], bed)
