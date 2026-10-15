#=
# Building a process from its parts

This example builds the four-step vacuum swing adsorption (VSA) process used
in [Getting started](getting_started.md) from its physical parameters,
following [Haghpanah et al. (2013)](https://doi.org/10.1021/ie302658y). A
process is described in layers, each of which can be edited, saved to a file
and reused on its own:

1. the gas and the feed;
2. the isotherm, given component by component, and the adsorbent;
3. the column;
4. the flowsheet: beds, and the devices connected to them;
5. the cycle: a list of steps;
6. the process, which puts them together with the time stepping.

Everything is in SI units. `Mocca.Units` has conversions, so that `0.1bar`
is `1e4` (Pa) and `-36.64kJ_per_mol` is `-36640.0` (J/mol).
=#

import Mocca
using Mocca.Units

# # 1. The gas and the feed
# Components are known by name. Their molar masses, heat capacities and
# viscosities default to standard values; here the molar mass of N₂ and the
# heat capacities follow Haghpanah et al. (2013). The mixture viscosity and
# density default to those of the feed. The first component is the one to
# capture.
gas = Mocca.Gas(["CO2", "N2"];
    molar_masses = (N2 = 28e-3,),
    heat_capacities = (CO2 = 697.5687, N2 = 1096.4),
    diffusivity = 1.6e-5)

feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = celsius(25), pressure = 1bar)

# # 2. The isotherm and the adsorbent
# A dual-site Langmuir isotherm. Each component lists its own parameters, so
# adding or removing a component does not mean editing parallel arrays.
# Capacities are per m³ of adsorbent here (`basis = :mass` takes mol/kg), and
# affinities are on a concentration basis [m³/mol]. `Mocca.langmuir` gives a
# single-site isotherm.
isotherm = Mocca.dual_site_langmuir(
    "CO2" => (qsb = 3489.44, b0 = 8.65e-7, ΔUb = -36.64121kJ_per_mol,
              qsd = 2872.35, d0 = 2.63e-8, ΔUd = -35.69066kJ_per_mol),
    "N2"  => (qsb = 6613.551, b0 = 2.5e-6, ΔUb = -15.8kJ_per_mol))

# The adsorbent has a density and heat capacity, and the particle properties
# that set the macropore diffusion of the linear driving force for mass
# transfer.
adsorbent = Mocca.Adsorbent("my zeolite 13X";
    isotherm = isotherm,
    density = 1130.0,              # kg/m³
    heat_capacity = 1070.0,        # J/(kg·K)
    particle_diameter = 2mm,
    particle_porosity = 0.35,
    tortuosity = 3.0)

# Check the isotherm before simulating: the loading of each component in
# equilibrium with the feed, in mol/kg, and the CO₂ isotherm at two
# temperatures.
q = Mocca.equilibrium_loading(adsorbent, gas; pressure = 1bar, temperature = 298.15, composition = (CO2 = 0.15, N2 = 0.85))
println("loading in equilibrium with the feed: CO₂ $(round(q[1], digits = 2)), N₂ $(round(q[2], digits = 3)) mol/kg")
pressures = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0] .* bar
table = Mocca.isotherm_table(adsorbent, gas, "CO2"; pressures = pressures, temperatures = [298.15, 348.15])
println(rpad("CO₂ at", 10), rpad("25 °C", 14), "75 °C")
for (P, row) in zip(pressures, eachrow(table))
    println(rpad("$(P / bar) bar", 10), join([rpad("$(round(x, digits = 2)) mol/kg", 14) for x in row]))
end

# The adsorbent can be saved and shared. This one is the zeolite 13X of
# Haghpanah et al. (2013), which Mocca's library also holds, as
# `Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")`.
path = Mocca.save_adsorbent(joinpath(mktempdir(), "my_13x.json"), adsorbent)
Mocca.load_adsorbent(path)

# # 3. The column
# The bed exchanges heat with the column wall, and the wall with the
# surroundings.
column = Mocca.Column(
    length = 1.0,                  # m
    radius = 14.45cm,
    outer_radius = 16.2cm,         # of the wall
    porosity = 0.37,               # bed void fraction
    axial_conductivity = 0.0903,   # W/(m·K)
    wall = Mocca.Wall(density = 7800.0, heat_capacity = 502.0, conductivity = 16.0),
    inner_heat_transfer = 8.6,     # W/(m²·K), bed to wall
    outer_heat_transfer = 2.5,     # W/(m²·K), wall to surroundings
    ambient_temperature = 298.15)

# # 4. The flowsheet
# A bed of 50 cells and three devices. Each device has a role, which sets its
# default boundary (the feed, the raffinate header, the vacuum line), tells
# the steps which device to open, and tells `performance` which streams are
# feed and product.
#
#     feed ──V_feed──┐                 ┌──V_product── raffinate
#                    bottom  Bed  top
#     extract ──V_vacuum──┘
#
# The bed has an energy balance for the wall by default (`thermal =
# Mocca.WithWall()`); `Mocca.Adiabatic()` and `Mocca.Isothermal()` are the
# alternatives (see the [breakthrough](breakthrough.md) example).
fs = Mocca.Flowsheet(gas, feed)
Mocca.add_bed!(fs, :Bed, adsorbent, column; ncells = 50)
Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))
Mocca.add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
fs

# # 5. The cycle
# Steps are written in process terms. Each opens the devices it needs, found
# by their roles, and keeps the others closed. Pressures carry over from step
# to step: the pressurisation starts where the evacuation left the bed, and
# the blowdown where the adsorption did. The feed velocity is interstitial.
# The other steps are `Purge`, `Equalise` and `Idle` (see the
# [two-bed VSA](two_bed_vsa.md) example), and new kinds of step can be added
# (see [Extending Mocca](@ref)).
cycle = Mocca.Cycle([
    Mocca.Pressurise(15.0),
    Mocca.Adsorb(15.0; velocity = 1.0),
    Mocca.Blowdown(30.0; to = 0.2bar),
    Mocca.Evacuate(40.0; to = 0.1bar),
])

# # 6. The process, and a check of the description
# By default a bed starts at the lowest pressure of its cycle; here it starts
# full of nitrogen at 1 bar. `validate` lists likely mistakes in the
# description that would not stop a simulation, such as a bed whose axial
# dispersion was set for a velocity far from that of its adsorption step.
process = Mocca.Process(fs, cycle; initial_state = Dict(:Bed => (pressure = 1bar,)))
issues = Mocca.validate(process)
println(isempty(issues) ? "no problems found" : join(sprint.(showerror, issues), "\n"))
process

# # 7. Simulate
# `simulate` runs the cycle a given number of times from the initial state.
# The outlet at the top of the bed over the first five cycles:
result = Mocca.simulate(process; cycles = 5)
bed_states = [s[:Bed] for s in result.states]
f_outlet = Mocca.plot_cell(bed_states, fs[:Bed], result.timesteps, 50)

# and the state along the bed at the end:
f_column = Mocca.plot_state(bed_states[end], fs[:Bed])

# The bed is still far from steady state after five cycles. The performance
# over the last of them, and at cyclic steady state:
perf_5 = Mocca.performance(result)
perf = Mocca.performance(Mocca.cyclic_steady_state(process; progress = true))
for (label, p) in ("after 5 cycles" => perf_5, "at steady state" => perf)
    println(rpad(label, 18), "purity $(round(100p.purity, digits = 1))%, recovery $(round(100p.recovery, digits = 1))%")
end

# # 8. Save the whole process
# The process file holds every layer, and reads back as the same process.
file = Mocca.write_process(joinpath(mktempdir(), "my_vsa.json"), process)
same = Mocca.read_process(file)
