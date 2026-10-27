#=
# Moving bed TSA for CO₂ capture from a waste-to-energy plant

This example follows [Mondino et al. (2022)](https://doi.org/10.1016/j.ijggc.2022.103685),
who design a moving bed temperature swing adsorption (MBTSA) process to capture
CO₂ from the flue gas of a 16.8 MWel waste-to-energy plant: 55.9 kg/s of gas
with 11% CO₂ at 30 °C, after drying.

In a moving bed the sorbent, here activated carbon beads, circulates through
sections that each do one job, instead of one column going through a cycle of
steps:

* **Adsorption:** flue gas flows up through the falling sorbent, which takes up
  the CO₂ (and some N₂); CO₂-free gas leaves the top.
* **Preheating:** a heat exchanger warms the sorbent with heat recovered from
  precooling. Most of the N₂ comes off and is recycled to the feed.
* **Desorption:** hot water at 187 °C heats the sorbent to about 184 °C; the CO₂
  comes off and is drawn from the bottom under a light vacuum (90 kPa).
* **Precooling and cooling:** the sorbent is cooled, first by the water loop
  that feeds preheating, then by cold water at 10 °C. A purge of CO₂-free gas
  flows up through both sections and is recycled to the feed.

```
sorbent ─► Adsorption ─► Preheating ─► Desorption ─► Precooling ─► Cooling ─┐
    ▲                                                                       │
    └───────────────────────────────────────────────────────────────────────┘
```

Each section is a moving bed on a [`Flowsheet`](@ref Mocca.Flowsheet), added
with [`add_moving_bed!`](@ref Mocca.add_moving_bed!). The sorbent passes from
one section to the next through [`connect_solids!`](@ref Mocca.connect_solids!),
the recovery water loop through
[`connect_heat_transfer_fluid!`](@ref Mocca.connect_heat_transfer_fluid!), and
the gas through flow devices, as for fixed beds. The process has no cycle of
steps: it is a [`SteadyProcess`](@ref Mocca.SteadyProcess), whose devices keep
their settings while it runs to steady state.

The sections, flows and temperatures are those of the paper's Table 6, and the
adsorbent is that of its Tables 2 and 5, which Mocca's library holds. The tube
area of the heat exchangers and the packing in the adsorption section come
from the SINTEF MBTSA MATLAB parameter set (`dcc3_template_AC.m`).
=#

import Mocca
using Mocca.Units
using CairoMakie

# # 1. The gas and the feed
# Dry flue gas, with the gas heat capacities of the paper. The light product,
# the CO₂-free gas that leaves the adsorption section, has a trace of CO₂.
gas = Mocca.Gas(["CO2", "N2"]; heat_capacities = (CO2 = 846.0, N2 = 1040.0), viscosity = 1.8e-5, diffusivity = 1.6e-5)
y_feed = [0.11, 0.89]
feed = Mocca.Feed(gas; composition = y_feed, temperature = celsius(30), pressure = 101.5e3,
    light_product = Mocca.light_product_composition(y_feed; trace = 1e-6))
feed_mass_flow = 55.9  # kg/s
molar_mass(y) = sum(y .* gas.molar_masses)
F_feed = feed_mass_flow / molar_mass(y_feed)

# # 2. The activated carbon
# Mondino et al. fit a Virial isotherm to CO₂ and N₂ measured on the activated
# carbon between 30 and 150 °C (their Table 5). Mocca's library holds it with
# the bead properties of their Table 2. Their macropores, 0.525 µm in radius,
# add Knudsen diffusion to the molecular diffusion that sets the rate of mass
# transfer:
ac = Mocca.load_adsorbent("activated_carbon_mondino_2022")

# The MATLAB parameter set has a competitive Langmuir isotherm for the same
# carbon, with loadings per kg. Compare the two at three points of the
# process: the feed end of the adsorption section, the CO₂ extraction point and
# the bottom of the cooling section.
ac_langmuir = Mocca.Adsorbent("activated carbon, MATLAB Langmuir";
    isotherm = Mocca.langmuir("CO2" => (qs = 0.4, b0 = 2.7e-7, ΔU = -27kJ_per_mol),
        "N2" => (qs = 0.02, b0 = 5.5e-7, ΔU = -18kJ_per_mol); basis = :mass),
    density = ac.density, heat_capacity = ac.heat_capacity, particle_diameter = ac.particle_diameter,
    particle_porosity = ac.particle_porosity, tortuosity = ac.tortuosity, macropore_radius = ac.macropore_radius)
conditions = [
    ("Feed, 30 °C, 11% CO₂", celsius(30), 101.5e3, y_feed),
    ("Extraction, 184 °C, 97% CO₂", celsius(184), 90e3, [0.97, 0.03]),
    ("Cooling, 18 °C, 1% CO₂", celsius(18), 101.3e3, [0.01, 0.99]),
]
println(rpad("", 30), "Virial (CO₂, N₂)      Langmuir (CO₂, N₂)   [mol/kg]")
for (label, T, P, y) in conditions
    qv = Mocca.equilibrium_loading(ac, gas; pressure = P, temperature = T, composition = y)
    ql = Mocca.equilibrium_loading(ac_langmuir, gas; pressure = P, temperature = T, composition = y)
    println(rpad(label, 30), rpad(string(round.(qv, digits = 3)), 22), round.(ql, sigdigits = 2))
end

# The Virial fit gives the loadings shown in the paper's Fig. 6: 0.35 mol/kg of
# CO₂ at the feed, falling to 0.07 at the extraction point. The Langmuir
# parameters give about a seventeenth of that CO₂ and almost no N₂. Their
# affinities `b0` (2.7e-7 and 5.5e-7 m³/mol) have the values of the Virial
# Henry constants at infinite temperature (2.6969e-7 and 5.5486e-7
# mol/(kg·kPa)), which suggests the Langmuir parameters were derived from
# Table 5 without converting between the two forms. The process below uses
# the Virial isotherm; section 9 repeats it with the Langmuir one.

# # 3. The sections
# Five sections of a stack, listed in the order the sorbent passes through
# them, with the heights, diameters and gas volume fractions of the paper's
# Table 6. The adsorption section also holds structured packing, 5% of its
# volume.
sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
heights = (0.9, 0.4, 0.8, 0.4, 0.6)
diameters = (12.8, 13.3, 13.3, 13.3, 13.3)
porosity = (0.7, 0.5, 0.5, 0.5, 0.5)
packing_fraction = (0.05, 0.0, 0.0, 0.0, 0.0)
sorbent_mass_flow = 650.0  # kg/s

# The axial dispersion and the effective axial conductivity of the gas follow
# the paper's correlations (Wakao and Funazkri 1978; Wakao et al. 1979, eqs.
# A.12 and A.21), `D_z = D_m/ε·(20 + 0.5·Sc·Re)` and `λ = k_g·(7 + 0.5·Pr·Re)`,
# at the feed velocity in the adsorption section and a slow flow of desorbed
# gas elsewhere.
area(i) = π * (diameters[i] / 2)^2
u_feed = F_feed * Mocca.GAS_CONSTANT * feed.temperature / feed.pressure / area(1)
velocity = (u_feed, 0.05, 0.05, 0.05, 0.05) ./ porosity
k_g, Pr, ρ_g, D_m, μ = 0.026, 0.7, 1.2, 1.6e-5, 1.8e-5
Re(i) = ρ_g * velocity[i] * porosity[i] * ac.particle_diameter / μ
Sc = μ / (ρ_g * D_m)
dispersion(i) = D_m / porosity[i] * (20 + 0.5 * Sc * Re(i))
conductivity(i) = k_g * (7 + 0.5 * Pr * Re(i))
columns = [Mocca.Column(length = heights[i], radius = diameters[i] / 2, porosity = porosity[i],
    axial_conductivity = conductivity(i)) for i in eachindex(sections)]

# Four sections have tube bundle heat exchangers, 229 m² of tube per m³ of
# section, with water in the tubes flowing up against the sorbent. The
# sorbent-side heat transfer coefficient follows the paper's correlation in the
# sorbent velocity (eq. 11). Desorption is heated with water at 187 °C and
# cooling cooled with water at 10 °C, both at 267.6 kg/s; preheating and
# precooling share a closed loop of 160.6 kg/s that carries heat from the hot
# sorbent leaving desorption to the sorbent entering it. Their inlet
# temperatures below are only starting values, since each takes the water
# leaving the other.
#
# The water-side heat transfer coefficient is not given in the paper. The
# demonstrator's water flow is laminar, with about 570 W/(m²·K) (Table 3); at
# full scale it is turbulent. With 10000 W/(m²·K) the heat loads of the
# preheating, desorption and cooling sections are within 1.5% of the paper's
# Table 7; with 3000 W/(m²·K), the Dittus–Boelter value for water at 0.3 m/s in
# the 6 mm bores, the heat recovered is 3% lower.
function heat_exchanger(i, fluid_mass_flow, fluid_inlet_temperature)
    solid_fraction = 1 - porosity[i] - packing_fraction[i]
    v_sorbent = sorbent_mass_flow / (ac.density * solid_fraction * area(i))
    return Mocca.tube_bundle(area_density = 229.0, sorbent_side_htc = Mocca.moving_bed_sorbent_htc(v_sorbent),
        fluid_side_htc = 10000.0, fluid_mass_flow = fluid_mass_flow, fluid_inlet_temperature = fluid_inlet_temperature)
end
heat_exchangers = (nothing,
    heat_exchanger(2, 160.6, feed.temperature),
    heat_exchanger(3, 267.6, celsius(187)),
    heat_exchanger(4, 160.6, feed.temperature),
    heat_exchanger(5, 267.6, celsius(10)))

# # 4. The flowsheet
# The flowsheet is built by a function of the adsorbent and the grid, so that
# section 8 can rebuild it on other grids. The gas flows through seven
# devices:
#
# * `V_feed` feeds the flue gas into the bottom of the adsorption section, at a
#   fixed rate.
# * `V_raffinate` lets the CO₂-free gas out of its top, through a valve to
#   atmospheric pressure, and `V_purge` takes part of it to the bottom of the
#   cooling section, at a fixed rate.
# * `V_link` passes the purge on from the top of cooling to the bottom of
#   precooling.
# * `V_recycle_preheat` and `V_recycle_precool` are blowers that return the gas
#   leaving the tops of preheating and precooling to the feed, holding those
#   tops at atmospheric pressure.
# * `V_product` draws the CO₂ from the bottom of desorption, through a valve to
#   90 kPa.
function mbtsa_flowsheet(adsorbent; cells_per_metre = 200)
    fs = Mocca.Flowsheet(gas, feed)
    for (i, name) in enumerate(sections)
        Mocca.add_moving_bed!(fs, name, adsorbent, columns[i];
            sorbent_mass_flow = sorbent_mass_flow, ncells = max(4, round(Int, cells_per_metre * heights[i])),
            packing_fraction = packing_fraction[i], heat_exchanger = heat_exchangers[i], dispersion = dispersion(i))
    end
    for (i, name) in enumerate(sections)
        Mocca.connect_solids!(fs, name => sections[mod1(i + 1, length(sections))])
    end
    Mocca.connect_heat_transfer_fluid!(fs, :Precooling => :Preheating)
    Mocca.connect_heat_transfer_fluid!(fs, :Preheating => :Precooling)

    Mocca.add_device!(fs, :V_feed, :feed; at = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_raffinate, :product; at = (:Adsorption, :top), boundary = (pressure = 101.325e3,))
    Mocca.add_device!(fs, :V_purge, :link; from = (:Adsorption, :top), to = (:Cooling, :bottom))
    Mocca.add_device!(fs, :V_link, :link; from = (:Cooling, :top), to = (:Precooling, :bottom))
    Mocca.add_device!(fs, :V_recycle_preheat, :link; from = (:Preheating, :top), to = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_recycle_precool, :link; from = (:Precooling, :top), to = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_product, :extract; at = (:Desorption, :bottom),
        boundary = (pressure = 90e3, temperature = celsius(187), composition = [0.97, 0.03]))
    return fs
end
fs = mbtsa_flowsheet(ac)

# # 5. The process
# The device settings. The purge is "approximately 15%" of the flue gas on a
# weight basis in the paper; here it is 15%, 8.4 kg/s. The cooling sorbent takes
# up about 0.3 mol/kg of N₂ from it, about 200 mol/s, and the sorbent's pores
# carry gas on to the adsorption section. What is left flows on through
# precooling, where it strips CO₂ from the hot sorbent and returns it to the
# feed. The purge has to be large enough to leave some: below about 7.9 kg/s
# none reaches the top of precooling, the sorbent goes back to the adsorption
# section with about 0.09 mol/kg of CO₂, and the capture rate falls to 87%.
# Above it the capture rate rises steeply, to 91.4% at 8.0 kg/s, and levels off
# at about 93% from 8.2 kg/s.
#
# The valves and blowers use the conductance of the ports they join (see
# [`default_conductance`](@ref Mocca.default_conductance)).
purge_mass_flow = 0.15 * feed_mass_flow  # kg/s
function mbtsa_settings(fs)
    C(d) = Mocca.default_conductance(fs, d)
    return Dict(
        :V_feed => Mocca.FixedMolarFlow(F_feed),
        :V_raffinate => Mocca.LinearValve(C(:V_raffinate)),
        :V_purge => Mocca.FixedMolarFlow(purge_mass_flow / molar_mass(feed.light_product)),
        :V_link => Mocca.LinearValve(C(:V_link)),
        :V_recycle_preheat => Mocca.BackPressureController(101.325e3, C(:V_recycle_preheat)),
        :V_recycle_precool => Mocca.BackPressureController(101.325e3, C(:V_recycle_precool)),
        :V_product => Mocca.LinearValve(C(:V_product)))
end

# The sections start at rough guesses of their steady temperatures, pressures
# and compositions, with the sorbent in equilibrium with the gas, which
# shortens the start-up. The process runs in spans of an hour, with steps
# growing to a minute, until the state changes little over a span; the
# sorbent goes round about fourteen times an hour.
initial = Dict(
    :Adsorption => (temperature = celsius(30), pressure = 101.5e3, composition = y_feed),
    :Preheating => (temperature = celsius(110), pressure = 101.325e3, composition = [0.5, 0.5]),
    :Desorption => (temperature = celsius(180), pressure = 90e3, composition = [0.95, 0.05]),
    :Precooling => (temperature = celsius(120), pressure = 101.325e3, composition = [0.5, 0.5]),
    :Cooling => (temperature = celsius(40), pressure = 101.325e3, composition = [0.02, 0.98]))
mbtsa_process(fs) = Mocca.SteadyProcess(fs, mbtsa_settings(fs); name = "MBTSA", duration = 3600.0,
    timestepping = Mocca.TimeStepping(max_dt = 60.0, first_dt = 0.05), heat_capacity_ratio = 1.4, initial_state = initial)
process = mbtsa_process(fs)

# # 6. Running to steady state
result = Mocca.steady_state(process; tol = 1e-3)

# The capture rate over the run shows the start-up settling within the first
# hour; the second hour checks that nothing changes any more:
capture(s) = Mocca.recovery(Mocca.stream_rates(s, :V_product), Mocca.stream_rates(s, :V_feed))
t = cumsum(result.timesteps) ./ 60
fig_startup = Figure(size = (700, 300))
ax = Axis(fig_startup[1, 1], xlabel = "Time [min]", ylabel = "CO₂ capture rate [%]")
lines!(ax, t, 100 .* capture.(result.states))
fig_startup

# The indicators of the paper's Table 7, at steady state. The paper's
# productivity is per tonne of sorbent in the process.
perf = Mocca.performance(result)
co2 = perf.extract[1] * gas.molar_masses[1]  # kg/s
rows = [
    ("CO₂ purity [%]", 100 * perf.purity, 97.2),
    ("CO₂ capture rate [%]", 100 * perf.recovery, 90.8),
    ("CO₂ captured [kg/s]", co2, 8.26),
    ("Specific heat duty [MJ/kg CO₂]", perf.specific_heat_duty / 1e6, 5.7),
    ("Productivity [kg CO₂/(t·h)]", co2 / perf.sorbent_inventory * 1e3 * 3600, 181),
    ("Cycle time [min]", perf.circulation_time / 60, 4.22),
    ("Preheating (recovered) [MW]", perf.heat_duty[:Preheating] / 1e6, 54.2),
    ("Desorption (external) [MW]", perf.heat_duty[:Desorption] / 1e6, 47.0),
    ("Precooling (recovered) [MW]", -perf.heat_duty[:Precooling] / 1e6, 54.2),
    ("Cooling (external) [MW]", -perf.heat_duty[:Cooling] / 1e6, 45.8),
]
println(rpad("", 34), "Mocca    Paper")
for (label, mocca, paper) in rows
    println(rpad(label, 34), rpad(round(mocca, digits = 2), 9), paper)
end

# The heat loads, the circulation and the productivity agree with the paper to
# within 2.5%. The capture rate is 2.3 points higher than the paper's, and the
# specific heat duty 3% lower. The water-side heat transfer coefficient is set
# so that the heat loads match (section 3); the purge is the paper's 15%
# (section 5). With it the sorbent brings about 0.03 mol/kg of CO₂ back to the
# adsorption section, against about 0.06 mol/kg in the paper's Fig. 6, which
# lets less CO₂ slip through. A purge just above the threshold of section 5
# gives the paper's capture rate, but the purity stays above 99%.
#
# The purity is 2.5 points higher than the paper's: the product holds about
# 0.3% N₂ against the paper's 2.8%. The models still differ in ways that may
# account for it: Mocca has one temperature for gas and sorbent, Darcy flow
# instead of Ergun's equation, and a single linear driving force for mass
# transfer, set by molecular and Knudsen diffusion in the macropores, without
# the paper's film and micropore resistances. Slower mass
# transfer, down to a tenth of the default, changes the capture rate and
# purity by less than 0.3 points, since the sorbent spends minutes in each
# section.
#
# The work: the feed compressor, the vacuum pump on the product and the fans
# and blowers that move the purge and the recycles up in pressure, and the
# work equivalent to the heating, the heat weighted by the Carnot factor at
# the temperatures the hot water gives it up at, against 25 °C:
work = [
    ("Compressor and vacuum pump [kW]", (perf.power - perf.blower_power) / 1e3),
    ("Fans and blowers [kW]", perf.blower_power / 1e3),
    ("Heating [MW]", perf.heating / 1e6),
    ("Equivalent work of the heating [MW]", perf.equivalent_work / 1e6),
    ("Equivalent work [MJ/kg CO₂]", 3.6 * perf.specific_equivalent_work / 1e3),
]
for (label, value) in work
    println(rpad(label, 38), round(value, sigdigits = 3))
end

# The electrical work is small next to the heat: about 160 kW for the feed
# compressor and the vacuum pump, and 8 kW for the fans and blowers, against
# 47 MW of heating. Weighted by its temperature, the heating is worth about
# 15 MW of work, a third of it, and the work equivalent to the whole process
# comes to about 1.8 MJ per kg of CO₂. [`add_measure!`](@ref Mocca.add_measure!)
# can weigh the heat differently, such as by the efficiency of the steam
# turbine the heat would otherwise go to.

# At steady state everything that enters leaves. The CO₂-free gas and the
# product together carry the feed:
state = result.states[end]
for (i, name) in enumerate(gas.names)
    n_out = perf.extract[i] + perf.product[i]
    println(rpad(name, 4), "in ", round(perf.feed[i], digits = 3), " mol/s, out ", round(n_out, digits = 3), " mol/s")
end

# The recycles return the N₂ released in preheating and the purge gas:
for d in (:V_recycle_preheat, :V_recycle_precool, :V_purge)
    n = Mocca.stream_rates(state, d)
    println(rpad(d, 18), round(sum(n), digits = 1), " mol/s, ", round(100 * n[1] / sum(n), digits = 1), "% CO₂")
end

# ## Profiles through the sections
# The layout of the paper's Fig. 6: from the bottom of the cooling section on
# the left to the top of the adsorption section on the right. The sorbent moves
# from right to left in each section.
fig_profiles = Mocca.plot_moving_bed_profiles(state, process.flowsheet)

# As in the paper, N₂ takes up a large share of the sorbent at the feed end of
# the adsorption section and comes off in preheating, which is why preheating
# gas is recycled. The CO₂ comes off in desorption. Precooling and cooling take
# the sorbent back down. The purge replaces the N₂ that the cold sorbent takes
# up, and the rest strips CO₂ from the hot sorbent in precooling, so the
# sorbent leaves precooling with about 0.03 mol/kg of CO₂, half the paper's;
# the gas grows rich in CO₂ towards the top of precooling. The temperatures at
# the ends of the sections, and those of the water, follow the paper's closely.
#
# In the adsorption section the CO₂ front is sharper than in the paper and sits
# near the top, with the sorbent below it close to equilibrium with the feed.

# # 7. Sensitivities
# How the capture rate and the work per tonne of CO₂ change with three
# operating settings: the sorbent flow, one variable that sets the flow of
# all five sections, the temperature of the hot water and the purge rate. A
# steady state is the fixed point of one long time step, so the adjoints that
# give the gradients of a cyclic process at its cyclic steady state give them
# here too, at the steady state found above (see
# [`performance_gradient`](@ref Mocca.performance_gradient)). The sorbent-side
# heat transfer coefficients were computed from the sorbent velocity when the
# flowsheet was built, and stay as they are.
Mocca.add_variable!(process, :sorbent_flow, ((:parameter, s, :SorbentMassFlow) for s in sections)...)
Mocca.add_variable!(process, :hot_water, (:parameter, :Desorption, :HXFluidInletTemperature))
Mocca.add_variable!(process, :purge, (:setting, :V_purge, :rate))
g = Mocca.performance_gradient(process; steady_state = result, measures = [:recovery, :specific_equivalent_work])

# As elasticities, the percentage change of each measure for a 1% increase of
# each setting (the hot water temperature in kelvin):
println(rpad("", 14), "capture   equivalent work")
for n in Mocca.variable_names(process)
    elasticity(m) = round(process[n] * g.gradients[m][n] / g.metrics[m], sigdigits = 3)
    println(rpad(n, 14), rpad(elasticity(:recovery), 10), elasticity(:specific_equivalent_work))
end

# More sorbent and hotter water both raise the capture rate, by about 0.4%
# for each 1% more, but at different costs. More sorbent raises the
# equivalent work by about 0.5% for each 1%, as there is more of it to heat.
# Water 1% hotter, about 4.6 K, raises it by over 4%: the heat is then given
# at a higher temperature, where it is worth more work. More purge raises the
# capture rate and lowers the work a little, as it strips more CO₂ from the
# hot sorbent in precooling (section 5). These are slopes at this design: the
# purge threshold of section 5 lies 6% below the purge used here, and beyond
# it the slopes change.

# # 8. Grid resolution
# What the sorbent and the water carry is advected with a limited second-order
# scheme, the van Albada limiter by default (see the `limiter` keyword of
# `add_moving_bed!`), which keeps fronts sharp. First-order upwinding,
# `limiter = Mocca.Upwind()`, smears the sorbent's heat in the cooling section
# as an axial conduction of about 10 W/(m·K) at 400 cells per metre, 50 times
# that of the gas; since this plant's capture rate hangs on how much purge gas
# gets through cooling, it then drifts by several points between 400 and 3200
# cells per metre. Each grid is run to steady state:
function steady_performance(adsorbent = ac; nonlinear_tolerance = 1e-6, kwarg...)
    p = mbtsa_process(mbtsa_flowsheet(adsorbent; kwarg...))
    return Mocca.performance(p; tol = 1e-3, nonlinear_tolerance = nonlinear_tolerance)
end
println("cells/m  purity  capture  MJ/kg")
for cpm in (25, 50, 100, 200, 400, 800)
    p = cpm == 200 ? perf : steady_performance(cells_per_metre = cpm)
    println(rpad(cpm, 9), rpad(round(100 * p.purity, digits = 1), 8), rpad(round(100 * p.recovery, digits = 1), 9),
        round(p.specific_heat_duty / 1e6, digits = 2))
end

# From 100 cells per metre the capture rate changes by about a point and the
# purity by less, and the default of 200 cells per metre is within 0.5 points
# of 800.

# # 9. With the Langmuir isotherm of the MATLAB parameter set
# The same plant with the Langmuir parameters captures little CO₂: the sorbent
# circulation was sized for a capacity seventeen times larger. The desorption
# gas is then almost pure CO₂, with the N₂ mole fraction at its lower bound,
# where the N₂ balance stops short of the default tolerance of 1e-6, so this
# run uses 1e-5.
perf_langmuir = steady_performance(ac_langmuir; nonlinear_tolerance = 1e-5)
println("Purity ", round(100 * perf_langmuir.purity, digits = 1), "%, capture rate ",
    round(100 * perf_langmuir.recovery, digits = 1), "%, specific heat duty ",
    round(perf_langmuir.specific_heat_duty / 1e6, digits = 0), " MJ/kg")

# # 10. Saving the process
# A steady process is saved like any other, with its moving beds, sorbent and
# water connections, and an `"operation"` with the device settings in place of
# a cycle:
path = Mocca.write_process(joinpath(mktempdir(), "mbtsa.json"), process)
Mocca.read_process(path)
