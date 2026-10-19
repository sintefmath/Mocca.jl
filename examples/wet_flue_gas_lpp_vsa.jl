#=
# Wet flue gas: single-bed VSA with light product pressurisation

This example follows section 3 of [Krishnamurthy et al.
(2014)](https://doi.org/10.1021/ie5024723), which captures CO₂ from a wet flue
gas of 15% CO₂, 82% N₂ and 3% H₂O at 25 °C with a single zeolite 13X bed. The
cycle has four steps:

* **Light product pressurisation (LPP):** nitrogen-rich light product enters
  the top of the evacuated bed and brings it up to the high pressure.
* **Adsorption:** feed enters the bottom at the high pressure; CO₂ and water
  are adsorbed and nitrogen leaves the top.
* **Blowdown:** the bed is taken down to an intermediate pressure through the
  top (co-current), which removes most of the remaining nitrogen.
* **Evacuation:** the bed is evacuated through the bottom (counter-current) to
  the low pressure, which gives the CO₂ product.

Water is adsorbed much more strongly than CO₂, so it stays in a short zone
near the feed end, but it takes up capacity there and makes the evacuation
cost more energy.

The model uses the paper's data: the isotherms (its Table 1, case 1), the
column, particle and gas properties (Table S2 of its Supporting Information),
the feed and the operating conditions. The particle size is not given there and
is taken from Haghpanah et al. (2013). Two simplifications: the LPP step lasts
a fixed 20 s, and uses pure nitrogen rather than the stored effluent of the
previous adsorption step.
=#

import Jutul
import Mocca

# # 1. Adsorbent
#
# The extended dual-site Langmuir isotherm on 13X for CO₂, N₂ and H₂O. At the
# feed conditions water fills the sorbent far more than CO₂ does:
isotherm = Mocca.wet_flue_gas_isotherm(:zeolite_13x)
C_feed = Mocca.WET_FLUE_GAS_FEED .* 101325.0 ./ (Mocca.GAS_CONSTANT * 298.15)
q_feed = Mocca.compute_equilibrium(isotherm, C_feed, 298.15)
for (name, q) in zip(Mocca.WET_FLUE_GAS_COMPONENTS, q_feed)
    println(rpad(name, 4), " loading in equilibrium with the feed: ", round(q, digits = 0), " mol/m³")
end

# # 2. Flowsheet
#
# `lpp_vsa_flowsheet` builds the bed and four devices, and returns the stages
# with the operating conditions of section 3.1 of the paper: adsorption 78.4 s,
# blowdown 37.9 s to 0.08 bar, evacuation 111.6 s to 0.03 bar, and a feed
# interstitial velocity of 0.45 m/s. The bed has 30 cells, as in the paper.
fs, stages = Mocca.lpp_vsa_flowsheet()
model = Mocca.setup_flowsheet_model(fs)
parameters = Mocca.setup_flowsheet_parameters(fs)
for st in stages
    println(rpad(st.name, 12), st.duration, " s: ", join(sort(collect(keys(st.settings))), ", "))
end

# # 3. Cycling to the paper's steady state
#
# The bed starts evacuated and full of nitrogen. The paper judges cyclic steady
# state by the CO₂ balance: the CO₂ fed and the CO₂ leaving must agree within
# 0.5% for five consecutive cycles, after at least 50 cycles (its eq. S19).
# `simulate_until_mass_balance` repeats the cycle until that holds.
#
# The time steps start at 0.01 s in every stage and double up to 1 s, which
# carries the solver through the abrupt start of a stage, such as the feed
# switching on at the start of adsorption.
state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.wet_flue_gas_initial_state(fs[:Bed], 0.03e5))
forces, timesteps = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0, first_dt = 0.01)
result = Mocca.simulate_until_mass_balance(model, state0, parameters, forces, timesteps;
    inflow = (:V_feed, :V_lpp), outflow = (:V_product, :V_vacuum), max_cycles = 300)
println("Cycles: $(result.cycles), CO₂ balance met: $(result.converged), last error: $(round(100 * result.history[end], digits = 2))%")

# The criterion only checks CO₂. Water is held far more strongly and loads the
# feed end of the bed over hundreds of cycles, slowly pushing CO₂ further along
# it, so the CO₂ balance drifts around 1% and may not close within 300 cycles.
# `Mocca.simulate_to_cyclic_steady_state` instead requires every variable,
# including the water loading, to settle.

# # 4. Process performance
#
# The metrics of the paper, from the streams through the devices over the last
# cycle: purity of the evacuation stream (with and without water), recovery of
# the CO₂ fed, vacuum pump energy for blowdown and evacuation per tonne of CO₂,
# and productivity per m³ of adsorbent per day.
states, ts = result.states, result.timesteps
M_CO2 = 44.01e-3
n_feed = Mocca.stream_totals(states, ts, :V_feed)
n_product = Mocca.stream_totals(states, ts, :V_vacuum)
energy = Mocca.vacuum_pump_energy(states, ts, :V_product) + Mocca.vacuum_pump_energy(states, ts, :V_vacuum)
cycle_time = sum(st -> st.duration, stages)
adsorbent_volume = sum(parameters[:Bed][:SolidVolume])
performance = (
    purity_wet = Mocca.purity(n_product),
    purity_dry = n_product[1] / (n_product[1] + n_product[2]),
    recovery = Mocca.recovery(n_product, n_feed),
    energy_kWh_per_t = Mocca.specific_energy_kwh_per_tonne(energy, n_product, M_CO2),
    productivity_t_per_m3_day = n_product[1] * M_CO2 / 1000 / adsorbent_volume / (cycle_time / 86400),
)
for (k, v) in pairs(performance)
    println(rpad(k, 26), round(v, sigdigits = 3))
end

# For the same conditions the paper reports 78.8% purity with water (95.2%
# without), 88.6% recovery, 188.8 kWh per tonne of CO₂ and 1.73 t of CO₂ per m³
# of 13X per day. Over 300 cycles this model gives about 83% purity with water
# (97.5% without), 79% recovery, 184 kWh/t and 1.55 t/m³/day, still drifting
# as water accumulates. Differences in the light product (nitrogen here, the
# stored effluent in the paper) and the numerical scheme contribute.

# # 5. Bed profiles
#
# Profiles at the end of each step of the last cycle, as in Figure 2 of the
# paper.
step_ends = cumsum([st.duration for st in stages])
f_profiles = Mocca.plot_profiles(states, ts, fs[:Bed], step_ends;
    labels = [st.name for st in stages], unit = :Bed, components = [1, 3])

# The streams through the devices over the cycle:
f_streams = Mocca.plot_flowsheet_streams(states, model, ts; stages = stages)
