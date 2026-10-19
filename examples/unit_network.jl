#=
# A flowsheet by hand

A [`Process`](@ref Mocca.Process) turns its device roles and steps into a
network of unit models, the device settings of each stage, and a Jutul
simulation case. This example builds the four-stage VSA cycle of Haghpanah et
al. (2013) directly from that layer, and checks that it gives the same results
as the process. A fixed bed is connected to three flow devices:

    feed ──V_feed──┐                 ┌──V_product── raffinate
                   bottom  Bed  top
    extract ──V_vacuum──┘

Each stage of the cycle is a set of device settings (valve open or closed,
flow controller rate, boundary pressure profile). Nothing inside the bed
changes between stages, so the same bed model can be wired into larger
flowsheets. This is the layer to work in when adding a kind of unit, device
law or step (see [Extending Mocca](@ref)).
=#

import Mocca
using Mocca.Units

# # 1. The process
# `four_stage_vsa_process` builds the flowsheet and its cycle of steps. The
# run uses steps of 1 s and starts from nitrogen at 1 bar.
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

#-
result = Mocca.simulate(process; cycles = 3, info_level = 0)
bed = process.flowsheet[:Bed]
bed_states = [s[:Bed] for s in result.states]

# # 2. The same flowsheet by hand
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

# # 3. Compare
# Both runs take the same time steps, so their states can be compared step by
# step:
y_process = [s[:y][1, ncells] for s in bed_states]
y_hand = [s[:Bed][:y][1, ncells] for s in states]
println("largest difference in the CO₂ fraction at the outlet: ", maximum(abs, y_process .- y_hand))

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
