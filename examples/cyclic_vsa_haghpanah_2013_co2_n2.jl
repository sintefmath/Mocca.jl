#=
# Cyclic Vacuum Swing Adsorption simulation
This example shows how to set up and run a cyclic vacuum swing adsorption
simulation as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y).

A two-component flue gas (CO₂ and N₂) is fed into a column of zeolite 13X,
which adsorbs CO₂ preferentially. The pressure in the column is then reduced
to desorb the CO₂ and collect it. The cycle has four steps:

* Pressurisation: the top of the column is closed and flue gas enters at the
  bottom, bringing the column up to the high pressure ``P_H``.
* Adsorption: flue gas is fed into the bottom at velocity ``v_{feed}``, and
  nitrogen-rich gas leaves the top at ``P_H``.
* Blowdown: the bottom is closed and the column is taken down to the
  intermediate pressure ``P_I`` through the top.
* Evacuation: the top is closed and the column is evacuated through the bottom
  to the low pressure ``P_L``, which gives the CO₂ product.

The process is described in four layers: the gas, the materials (adsorbent
and column), the flowsheet (a bed and the devices connected to it) and the
cycle (a list of steps). Adsorption on zeolite 13X follows a dual-site
Langmuir isotherm.
=#

import Mocca
using Mocca.Units

# # 1. Gas and feed
# CO₂ and N₂ with the properties used by Haghpanah et al. Components are known
# to Mocca by name; the properties given here override its defaults. The
# first component is the one to capture.
gas = Mocca.Gas(["CO2", "N2"];
    molar_masses = [44.01e-3, 28e-3],          # kg/mol
    heat_capacities = [697.5687, 1096.4],      # J/(kg·K)
    viscosity = 1.72e-5,                       # Pa·s
    density = 1.22638310956,                   # kg/m³
    diffusivity = 1.6e-5)                      # m²/s

feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)

# # 2. Adsorbent and column
# The zeolite 13X of the paper is in Mocca's library of adsorbents; the
# [custom setup](custom_setup_cyclic_vsa.md) example builds one from its
# parameters.
adsorbent = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")

column = Mocca.Column(length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    axial_conductivity = 0.0903,
    wall = Mocca.Wall(density = 7800.0, heat_capacity = 502.0, conductivity = 16.0),
    inner_heat_transfer = 8.6, outer_heat_transfer = 2.5, ambient_temperature = 298.15)

# # 3. Flowsheet
# A bed of 200 cells with three devices. Each device has a role, which sets
# its default boundary (the feed, the raffinate header, the vacuum line) and
# tells the steps which device to open:
#
#     feed ──V_feed──┐                 ┌──V_product── raffinate
#                    bottom  Bed  top
#     extract ──V_vacuum──┘
fs = Mocca.Flowsheet(gas, feed)
Mocca.add_bed!(fs, :Bed, adsorbent, column; ncells = 200)
Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))
Mocca.add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))

# # 4. Cycle
# Pressures carry over from step to step: pressurisation starts where the
# evacuation left the bed, and blowdown where adsorption did. The feed
# velocity is interstitial.
cycle = Mocca.Cycle([
    Mocca.Pressurise(15.0),
    Mocca.Adsorb(15.0; velocity = 1.0),
    Mocca.Blowdown(30.0; to = 0.2bar),
    Mocca.Evacuate(40.0; to = 0.1bar),
])

# # 5. Process
# The column starts full of nitrogen at 1 bar. By default a bed starts at the
# lowest pressure of its cycle.
process = Mocca.Process(fs, cycle; initial_state = Dict(:Bed => (pressure = 1bar,)))

# The same process is ready-made as `Mocca.four_stage_vsa_process()`.
#
# # 6. Simulate
# We run 3 cycles for demonstration purposes; `Mocca.cyclic_steady_state`
# finds the steady state.
result = Mocca.simulate(process; cycles = 3, info_level = 0)

# # 7. Visualisation
# The primary variables at the outlet (top) of the bed through time:
bed = fs[:Bed]
bed_states = [s[:Bed] for s in result.states]
f_outlet = Mocca.plot_cell(bed_states, bed, result.timesteps, 200)

# and along the column at the end of the simulation:
f_column = Mocca.plot_state(bed_states[end], bed)

# The streams through the devices:
f_streams = Mocca.plot_flowsheet_streams(result.states, process.model, result.timesteps; stages = Mocca.stages(process))

# The performance over the last cycle. It is still far from steady state
# after three cycles.
perf = Mocca.performance(result)
println("purity $(round(100perf.purity, digits = 1))%, recovery $(round(100perf.recovery, digits = 1))%")
