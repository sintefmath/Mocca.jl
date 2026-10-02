# # Direct Column Breakthrough simulation

# This example shows how to set up and run a direct column breakthrough (DCB)
# simulation as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y).
# A two-component flue gas (CO₂ and N₂) is fed into a column of zeolite 13X
# initially filled with N₂, until CO₂ breaks through at the outlet.
#
# Adsorption on zeolite 13X follows a dual-site Langmuir isotherm. The feed
# enters the bottom of the column at a fixed velocity, and the top is open at
# the feed pressure. There is no heat transfer between the column and its
# wall.

import Mocca
using Mocca.Units

# # Gas, adsorbent and column
# The gas and adsorbent are those of the paper; see the
# [cyclic VSA](cyclic_vsa_haghpanah_2013_co2_n2.md) example. For the DCB the
# heat transfer coefficients to the wall are zero.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)
adsorbent = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")
column = Mocca.Column(length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    axial_conductivity = 0.0903, inner_heat_transfer = 0.0, outer_heat_transfer = 0.0)

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
# column.
ncells = 200
fs = Mocca.Flowsheet(gas, feed)
Mocca.add_bed!(fs, :Bed, adsorbent, column; ncells = ncells)
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
result = Mocca.simulate(process; timestep_selector_cfg = (Bed = (y = 0.01, Temperature = 10.0, Pressure = 10.0),))

# # Visualisation
# The primary variables at the outlet through time:
bed = fs[:Bed]
bed_states = [s[:Bed] for s in result.states]
f_outlet = Mocca.plot_cell(bed_states, bed, result.timesteps, ncells)

# and along the column at the end of the simulation:
f_column = Mocca.plot_state(bed_states[end], bed)
