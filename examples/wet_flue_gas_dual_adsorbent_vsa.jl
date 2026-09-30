#=
# Wet flue gas: dual-adsorbent, two-bed VSA

This example sets up the process proposed in section 4 of [Krishnamurthy et
al. (2014)](https://doi.org/10.1021/ie5024723) for CO₂ capture from a wet flue
gas of 15% CO₂, 82% N₂ and 3% H₂O. A silica gel bed takes up the water, and
the dried gas passes to a zeolite 13X bed that concentrates the CO₂. Keeping
the adsorbents in separate beds, rather than layering them in one, lets each
be evacuated on its own, so the CO₂ product comes out practically dry.

    feed ──V_feed──┐                      ┌──V_link──┐                   ┌──V_product── raffinate
                   bottom  SilicaGel  top ┘          bottom  Zeolite  top┤
    H2O ◄─V_waste──┘                                 │                   └──V_lpp────── light product
                                          CO2 ◄─V_vacuum┘

Both beds run a four-step cycle, coupled only during adsorption:

| Step           | Silica gel bed                       | 13X bed                                   |
|----------------|--------------------------------------|-------------------------------------------|
| pressurisation | with feed, from the bottom           | with light product, from the top          |
| adsorption     | feed in; effluent to the 13X bed     | silica gel effluent in; nitrogen out      |
| blowdown       | out of the bottom (counter-current)  | out of the top (co-current)               |
| evacuation     | out of the bottom: the water         | out of the bottom: the CO₂ product        |

The blowdown and evacuation times and pressures differ between the beds. The
silica gel bed finishes first and stays closed until the 13X bed completes
its cycle.

The model uses the paper's data: the isotherms (its Table 1, case 1, for both
adsorbents), the column, particle and gas properties (Table S2 of its
Supporting Information), the feed and the operating conditions for minimum
energy (section 4.3). The particle size is not given there and is taken from
Haghpanah et al. (2013). Two simplifications: the 13X light product pressurisation lasts a fixed 20 s,
like the silica gel pressurisation, and uses pure nitrogen rather than the
stored effluent of the previous adsorption step.
=#

import Jutul
import Mocca

# # 1. Adsorbents
#
# Silica gel adsorbs little CO₂ and N₂ but a great deal of water; 13X adsorbs
# both CO₂ and water strongly. Loadings in equilibrium with the feed:
C_feed = Mocca.WET_FLUE_GAS_FEED .* 101325.0 ./ (Mocca.GAS_CONSTANT * 298.15)
for adsorbent in (:silica_gel, :zeolite_13x)
    q = Mocca.compute_equilibrium(Mocca.wet_flue_gas_isotherm(adsorbent), C_feed, 298.15)
    println(rpad(adsorbent, 12), join(["$n $(round(qi, digits = 0))" for (n, qi) in zip(Mocca.WET_FLUE_GAS_COMPONENTS, q)], ", "), " mol/m³")
end

# # 2. Flowsheet
#
# `dual_adsorbent_vsa_flowsheet` builds the two beds and six devices. Its
# defaults are the operating conditions for minimum energy in section 4.3 of
# the paper:
#
# | | Silica gel bed | 13X bed |
# |---|---|---|
# | length | 0.41 m | 1 m |
# | pressurisation | 20 s | 20 s |
# | adsorption | 46.2 s | 46.2 s |
# | blowdown | 42.21 s to 0.48 bar | 56.3 s to 0.07 bar |
# | evacuation | 46.02 s to 0.3 bar | 101.2 s to 0.03 bar |
#
# with a feed interstitial velocity of 0.7 m/s. Each bed has 30 cells, as in
# the paper. The per-bed steps are combined into one sequence of stages, split
# wherever either bed changes step:
fs, stages = Mocca.dual_adsorbent_vsa_flowsheet()
model = Mocca.setup_flowsheet_model(fs)
parameters = Mocca.setup_flowsheet_parameters(fs)
for st in stages
    println(rpad(st.name, 42), round(st.duration, digits = 2), " s: ", join(sort(collect(keys(st.settings))), ", "))
end

# # 3. Cycling to the paper's steady state
#
# As in the paper, the silica gel bed starts at its low pressure and the 13X
# bed at the high pressure, both full of nitrogen. The paper judges cyclic
# steady state by the CO₂ balance: the CO₂ fed and the CO₂ leaving must agree
# within 0.5% for five consecutive cycles, after at least 50 cycles (its eq.
# S19). `simulate_until_mass_balance` repeats the cycle until that holds.
#
# The time steps start at 0.01 s in every stage and double up to 1 s, which
# carries the solver through the abrupt start of a stage, such as the feed
# switching on at the start of adsorption.
state0 = Mocca.setup_flowsheet_state(fs;
    SilicaGel = Mocca.wet_flue_gas_initial_state(fs[:SilicaGel], 0.30e5),
    Zeolite = Mocca.wet_flue_gas_initial_state(fs[:Zeolite], 101325.0))
forces, timesteps = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0, first_dt = 0.01)
result = Mocca.simulate_until_mass_balance(model, state0, parameters, forces, timesteps;
    inflow = (:V_feed, :V_lpp), outflow = (:V_waste, :V_product, :V_vacuum), max_cycles = 300)
println("Cycles: $(result.cycles), CO₂ balance met: $(result.converged), last error: $(round(100 * result.history[end], digits = 2))%")

# The criterion only checks CO₂, and the water is far from steady when it is
# met. The silica gel keeps about 98% of the water fed each cycle and removes
# little of it in blowdown and evacuation, so it fills over the cycles. In this
# model the water starts to break through to the 13X bed after about 100
# cycles, close to when the CO₂ balance closes, and rises from about 20 ppm at
# cycle 100 to over 1% by cycle 300. `Mocca.simulate_to_cyclic_steady_state`
# requires every variable, including the water loading, to settle; for this
# process it reaches a state in which the silica gel is saturated and much of
# the water passes into the 13X bed.

# # 4. Process performance
#
# The metrics of the paper (its equations 9–14), from the streams through the
# devices over the last cycle. The CO₂ product is the 13X evacuation stream;
# recovery is relative to the CO₂ fed to the silica gel bed; the energy is the
# vacuum pump work for blowdown and evacuation of both beds per tonne of CO₂;
# productivity is per m³ of 13X per day.
states, ts = result.states, result.timesteps
M_CO2 = 44.01e-3
n_feed = Mocca.stream_totals(states, ts, :V_feed)
n_product = Mocca.stream_totals(states, ts, :V_vacuum)
n_to_zeolite = Mocca.stream_totals(states, ts, :V_link)
energy = sum(Mocca.vacuum_pump_energy(states, ts, d) for d in (:V_waste, :V_product, :V_vacuum))
cycle_time = sum(st -> st.duration, stages)
zeolite_volume = sum(parameters[:Zeolite][:SolidVolume])
performance = (
    purity = Mocca.purity(n_product),
    recovery = Mocca.recovery(n_product, n_feed),
    energy_kWh_per_t = Mocca.specific_energy_kwh_per_tonne(energy, n_product, M_CO2),
    productivity_t_per_m3_day = n_product[1] * M_CO2 / 1000 / zeolite_volume / (cycle_time / 86400),
    water_to_13X_fraction_of_feed = n_to_zeolite[3] / n_feed[3],
    water_to_13X_max_ppm = 1e6 * maximum(s[:V_link][:InletComposition][3, 1] for s in states if s[:V_link][:MolarFlow][1] > 0),
)
for (k, v) in pairs(performance)
    println(rpad(k, 30), round(v, sigdigits = 3))
end

# For these conditions the paper reports 95% purity, 90% recovery, 177 kWh per
# tonne of CO₂ and 1.82 t of CO₂ per m³ of 13X per day, with 20 ppm of water in
# the gas entering the 13X bed. When the CO₂ balance closes, this model gives
# about 95.6% purity, 86% recovery, 190 kWh/t and 1.75 t/m³/day, with a few
# hundred ppm of water entering the 13X bed.

# # 5. Bed profiles
#
# Profiles at the end of each bed's steps in the last cycle, as in Figures 10
# and 11 of the paper.
t_press, t_ads = 20.0, 46.20
sg_ends = cumsum([t_press, t_ads, 42.21, 46.02])
x_ends = cumsum([t_press, t_ads, 56.30, 101.20])
step_names = ["pressurisation", "adsorption", "blowdown", "evacuation"]
f_silica_gel = Mocca.plot_profiles(states, ts, fs[:SilicaGel], sg_ends;
    labels = step_names, unit = :SilicaGel, components = [1, 3])
f_zeolite = Mocca.plot_profiles(states, ts, fs[:Zeolite], x_ends;
    labels = step_names, unit = :Zeolite, components = [1, 3])

# The streams through the devices over the cycle:
f_streams = Mocca.plot_flowsheet_streams(states, model, ts; stages = stages)
