# # Breakthrough curve

# This example simulates the direct column breakthrough (DCB) of
# [Haghpanah et al. (2013)](http://dx.doi.org/10.1021/ie302658y), and how the
# heat released by adsorption shapes it. A two-component flue gas (CO₂ and
# N₂) is fed into a column of zeolite 13X initially filled with N₂, until CO₂
# breaks through at the outlet. A breakthrough is a single step rather than a
# cycle, so it is run once with `simulate`.

import Mocca
using Mocca.Units

# # Gas, adsorbent and column
# The gas, adsorbent and column of the paper; the
# [process from its parts](process_from_parts.md) example sets them up from
# their parameters.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)
adsorbent = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")
column = Mocca.haghpanah_2013_column()

#-
# The gas velocity follows Darcy's law,
#
# ```math
# v = -\frac{k}{\mu}(\nabla P)
# ```
#
# with the permeability `k` of a packed bed of particles of diameter ``d_p``
# and bed porosity ``\varepsilon`` from the Kozeny–Carman relation,
#
# ```math
# k = \frac{4}{150}\frac{\varepsilon^3}{(1-\varepsilon)^2} \left(\frac{d_p}{2}\right)^2.
# ```
#
# `add_bed!` computes it, and the axial dispersion, from the adsorbent and
# column. In the DCB of the paper no heat leaves the column, so the bed is
# adiabatic.
ncells = 200
fs = Mocca.Flowsheet(gas, feed)
Mocca.add_bed!(fs, :Bed, adsorbent, column; ncells = ncells, thermal = Mocca.Adiabatic())
Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))

# # The breakthrough as a one-step cycle
# A single adsorption step of 5000 s at an interstitial velocity of 1 m/s.
# The column starts at the feed pressure, the lowest pressure of this cycle,
# full of nitrogen.
cycle = Mocca.Cycle([Mocca.Adsorb(5000.0; velocity = 1.0)])
process = Mocca.Process(fs, cycle; timestepping = Mocca.TimeStepping(max_dt = 5000.0, first_dt = 5000.0))

# # Run the simulation
# One report step of 5000 s, which the solver divides according to target
# changes of the state variables of the bed per step.
selector = (Bed = (y = 0.01, Temperature = 10.0, Pressure = 10.0),)
result = Mocca.simulate(process; timestep_selector_cfg = selector)

# # Visualisation
# The primary variables at the outlet through time. CO₂ reaches the outlet in
# two steps. The heat released by adsorption travels with the first front, so
# the bed behind it is about 60 K hotter and holds less CO₂; the rest of the
# CO₂ breaks through as the feed cools the bed again:
bed = fs[:Bed]
bed_states = [s[:Bed] for s in result.states]
f_outlet = Mocca.plot_cell(bed_states, bed, result.timesteps, :top)

# and along the column at the end of the simulation:
f_column = Mocca.plot_state(bed_states[end], bed)

# # Heat effects
# Each bed has a thermal model: `Mocca.Adiabatic()`, `Mocca.WithWall()` (the
# default), which adds an energy balance for the column wall and its heat
# exchange with the gas and the surroundings, or `Mocca.Isothermal(T)`, with
# no energy balance. The same breakthrough with each:
function breakthrough(thermal)
    fs = Mocca.Flowsheet(gas, feed)
    Mocca.add_bed!(fs, :Bed, adsorbent, column; ncells = ncells, thermal = thermal)
    Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))
    p = Mocca.Process(fs, cycle; timestepping = process.timestepping)
    r = Mocca.simulate(p; timestep_selector_cfg = selector)
    return [s[:Bed] for s in r.states], r.timesteps
end

runs = [
    "Adiabatic" => (bed_states, result.timesteps),
    "With wall" => breakthrough(Mocca.WithWall()),
    "Isothermal" => breakthrough(Mocca.Isothermal(298.15)),
]

# The time at which the outlet reaches half the CO₂ fraction of the feed:
for (label, (states, ts)) in runs
    k = findfirst(s -> s[:y][1, end] > 0.075, states)
    println(rpad(label, 12), round(cumsum(ts)[k], digits = 0), " s")
end

# The isothermal bed stays at the feed temperature and holds more CO₂, so CO₂
# breaks through about twice as late, in a single sharp front. Heat lost
# through the wall cools the bed a little behind the first front, which then
# takes up more CO₂ than the adiabatic bed.
f_comparison = Mocca.plot_cell_comparison(runs, bed, :top)
