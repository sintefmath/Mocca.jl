#=
# Direct air capture by steam-assisted TVSA

This example follows [Stampi-Bombelli et al.
(2020)](https://doi.org/10.1007/s10450-020-00249-w), who capture CO₂ from
ambient air on an amine-functionalised nanofibrillated cellulose (APDES-NFC)
by temperature–vacuum swing adsorption with a steam purge (S-TVSA). The air
holds 400 ppm of CO₂ and is at 20 °C and 50% relative humidity. The bed is a
pancake, 1 cm long and 8 cm across, inside a heating jacket. The cycle has
four steps:

* **Adsorption:** air flows through the bed at 1 bar for almost four hours,
  until the bed is close to saturated with CO₂. The sorbent also takes up
  water, which raises how much CO₂ it holds.
* **Blowdown:** a vacuum pump at the top takes the bed down to 50 mbar in
  30 s.
* **Heating:** the jacket is heated to 95 °C while the pump holds the bed at
  50 mbar, which drives the remaining air out.
* **Desorption:** steam at 95 °C and 50 mbar flows up through the bed for over
  eight hours, while the pump draws CO₂ and water out of the top. The water is
  condensed before the pump.

Three things in Mocca make this process: the step's `jacket_temperature`,
which heats the wall during a step, the paper's isotherm for CO₂ and water
([`toth_gab`](@ref)) and mass transfer coefficients fitted to breakthrough
curves ([`lumped_ldf`](@ref)).

The model uses the paper's data: the isotherms of its Table 2, the column,
sorbent and transfer coefficients of its Table 1 and the operating conditions
of its Table 5. Where it differs:

* a 30 s pressurisation with air through the feed end comes before
  adsorption. The paper opens the bed to 1 bar at the start of adsorption,
  which is too abrupt for the solver;
* the CO₂ isotherm depends on the water loading in equilibrium with the gas,
  where the paper does not say which water loading it uses;
* the paper gives no permeability, only a pressure drop below 1 mbar; the
  Blake–Kozeny permeability of its particles gives 0.05 mbar;
* the adsorbed phase takes the heat capacity of its components weighted by
  the composition of the gas, not of the adsorbed phase.
=#

import Mocca
using Mocca.Units
using CairoMakie

# # 1. Gas and adsorbent
#
# Air is CO₂, water and nitrogen, as in the paper. The sorbent comes from
# Mocca's library. At the feed conditions it holds about 2.4 mol/kg of CO₂, and
# without the water only 2.3: water raises the affinity of the Toth isotherm
# for CO₂ and its capacity. At 95 °C in steam with 5% CO₂ it holds a few
# percent of that, which is the swing the cycle works with.
gas = Mocca.stampi_bombelli_2020_gas()
i_N2, i_H2O = findfirst(==("N2"), gas.names), findfirst(==("H2O"), gas.names)
adsorbent = Mocca.load_adsorbent("apdes_nfc_stampi_bombelli_2020")
air = (CO2 = 0.0004, H2O = 0.0115, N2 = 0.9881)
for (label, P, T, y) in [("humid air, 20 °C", 1bar, 293.15, air), ("dry air, 20 °C", 1bar, 293.15, (CO2 = 0.0004, N2 = 0.9996)),
        ("steam with 5% CO₂, 95 °C", 50mbar, 368.15, (CO2 = 0.05, H2O = 0.95))]
    q = Mocca.equilibrium_loading(adsorbent, gas; pressure = P, temperature = T, composition = y)
    println(rpad(label, 26), "CO₂ ", round(q[1], digits = 3), " mol/kg, H₂O ", round(q[i_H2O], digits = 3), " mol/kg")
end

# # 2. Process
#
# [`dac_tvsa_process`](@ref) builds the bed in its jacket, the devices and the
# cycle at the base case of the paper: air at 50 cm³/s, desorption at 50 mbar
# and 95 °C, steam at 25 cm³/s, and the step times of its Table 5. The bed has
# 10 cells; 5 or 40 change the CO₂ produced by less than 0.2%, since the slow
# mass transfer spreads the fronts over the whole bed.
process = Mocca.dac_tvsa_process()

# The steps of its cycle, written out. The jacket is at ambient temperature
# unless a step sets `jacket_temperature`. Heating is an evacuation that holds
# the bed at 50 mbar through the vent, and desorption a purge with steam
# through the light inlet at the bottom, out through the extract at the top:
#
# ```julia
# Mocca.Cycle([
#     Mocca.Pressurise(30.0),
#     Mocca.Adsorb(13772.0; flow_rate = 50e-6),
#     Mocca.Blowdown(30.0; to = 50mbar, device = :V_vent),
#     Mocca.Evacuate(704.0; to = 50mbar, through = :top, device = :V_vent, jacket_temperature = 368.15, name = "heating"),
#     Mocca.Purge(30000.0; flow_rate = 25e-6, jacket_temperature = 368.15, name = "desorption"),
# ])
# ```

# # 3. Cyclic steady state
#
# Each adsorption nearly saturates the bed and each desorption nearly empties
# it, so the bed forgets where it started within a few cycles: three plain
# cycles bring the change over a cycle to 10⁻⁶, and Newton's method takes it
# below 10⁻⁹ in one or two iterations.
result = Mocca.cyclic_steady_state(process; progress = true)
println("Converged: $(result.converged)")

# # 4. Performance
#
# The paper's measures, from the streams through the devices over the last
# cycle:
#
# * the CO₂ production rate, in kg/h, and its purity without water, which is
#   condensed out before the pump;
# * the thermal energy: the heat from the jacket ([`jacket_heat`](@ref), which
#   [`performance`](@ref) reports as `heating`) and the heat to raise the
#   steam from water at ambient temperature, its eq. (19);
# * the electrical energy: the fan that blows air through the bed, and the
#   vacuum pump, which compresses the dry gas isothermally from the bed to
#   1 bar, with the efficiencies of its Table 4 (eqs. 20 to 22).
#
# Both energies are per kg of CO₂ produced.
function dac_performance(result; T_amb = 293.15, P_H = 1bar, η_fan = 0.5)
    p = result.process
    perf = Mocca.performance(result)
    states, ts = result.states, result.timesteps
    M = p.flowsheet.gas.molar_masses
    R = 8.314462618
    n = perf.extract
    m_CO2 = n[1] * M[1]
    ## Steam, raised from water at T_amb to the boiling point at P_L and
    ## superheated to its temperature
    steam = Mocca.boundary_condition(p.flowsheet, :V_steam, :inlet)
    P_L, T_s = steam.pressure, steam.temperature
    T_b = 1730.63 / (8.07131 - log10(P_L / 133.322368)) - 233.426 + 273.15
    ΔH_v = 2.501e6 - 2370.0 * (T_b - 273.15)
    m_steam = Mocca.stream_totals(states, ts, :V_steam)[i_H2O] * M[i_H2O]
    Q_steam = m_steam * (4186.0 * (T_b - T_amb) + ΔH_v + 1864.0 * (T_s - T_b))
    Q = max(0.0, perf.heating) + Q_steam
    ## Fan: the air flow times the pressure drop over the bed
    W_fan = 0.0
    for (s, dt) in zip(states, ts)
        F = s[:V_feed][:MolarFlow][1]
        F > 0 || continue
        P = s[:Bed][:Pressure]
        V = F * R * s[:V_feed][:InletTemperature][1] / P[1]
        W_fan += V * (P[1] - P[end]) / η_fan * dt
    end
    ## Vacuum pump, on the dry gas, with an efficiency from 0.3 at 10 mbar to
    ## 0.7 at 100 mbar
    η_pump = clamp(0.3 + 0.4 * (P_L - 10mbar) / (100mbar - 10mbar), 0.3, 0.7)
    W_pump = 0.0
    for d in (:V_vent, :V_extract), (s, dt) in zip(states, ts)
        F, P_out = s[d][:MolarFlow][1], s[d][:OutletPressure][1]
        (F > 0 && P_out < P_H) || continue
        F_dry = F * (1 - s[d][:InletComposition][i_H2O, 1])
        W_pump += F_dry * R * T_amb * log(P_H / P_out) / η_pump * dt
    end
    return (production = m_CO2 / perf.cycle_time * 3600, purity_dry = n[1] / (n[1] + n[i_N2]),
        recovery = perf.recovery, thermal = Q / m_CO2 / 1e6, electrical = (W_fan + W_pump) / m_CO2 / 1e6,
        jacket = perf.heating / m_CO2 / 1e6, steam = Q_steam / m_CO2 / 1e6)
end

base = dac_performance(result)
println("CO₂ production       ", round(base.production * 1e5, sigdigits = 3), "×10⁻⁵ kg/h")
println("CO₂ purity, dry      ", round(100 * base.purity_dry, digits = 2), "%")
println("CO₂ recovered        ", round(100 * base.recovery, digits = 1), "% of the CO₂ in the air")
println("thermal energy       ", round(base.thermal, sigdigits = 3), " MJ/kg (jacket ", round(base.jacket, sigdigits = 3),
    ", steam ", round(base.steam, sigdigits = 3), ")")
println("electrical energy    ", round(base.electrical, sigdigits = 3), " MJ/kg")

# The paper produces 1.2 to 2.6×10⁻⁵ kg/h over its desorption times, at a dry
# purity above 99.6%, and this model gives 2.2×10⁻⁵ kg/h at 99.96%. About half
# the CO₂ in the air is captured: the adsorption runs until the bed is close to
# saturated, and the rest passes through. The electrical energy, 0.41 MJ/kg,
# is above the paper's 0.25 to 0.3 MJ/kg at 50 mbar; most of it is the vacuum
# pump during blowdown and desorption.
#
# The thermal energy is mostly the steam. Raising 25 cm³/s of steam at 50 mbar
# and 95 °C for 30000 s takes 58 kJ per cycle, about 210 MJ per kg of CO₂,
# where the paper's Figure 6 shows 35 to 100 MJ/kg for this steam flow. Its
# steam flow seems to be measured at other conditions, which the paper does not
# give. The heat from the jacket, which warms the bed and wall and desorbs CO₂
# and water, is 6 MJ/kg.

# # 5. Through the cycle at the top of the bed
#
# As Figure 4 of the paper: the state of the last cell, at the top of the bed,
# through the last cycle. CO₂ is first taken up near the feed end, and reaches
# the top after about 40 minutes. Its loading there then rises to 2.15 mol/kg,
# against 1.8 in the paper, and is still rising when adsorption ends. The
# paper lowered the capacities of its isotherms to match its breakthrough
# experiments, by an amount it does not give, which likely explains why this
# model holds more CO₂, and produces about a quarter more over the cycle.
# Blowdown takes out part of the air and heating the rest, while the bed warms
# to 50 °C. Steam then carries the CO₂ out, and the bed passes 90 °C within
# about two hours and approaches 95 °C over the rest of the desorption. The water loading follows the humidity of the gas: up to
# 2.8 mol/kg at 50% humidity in adsorption, and down to 0.2 mol/kg in the hot
# steam, at 6% humidity, as in the paper.
function exit_figure(result)
    p = result.process
    ρ = Mocca.load_adsorbent("apdes_nfc_stampi_bombelli_2020").density
    bed = [s[:Bed] for s in result.states]
    t = cumsum(result.timesteps)
    t .-= t[1] - result.timesteps[1]
    T = [s[:Temperature][end] for s in bed]
    series = [
        ("CO₂ loading [mol/kg]", [s[:AdsorbedConcentration][1, end] / ρ for s in bed]),
        ("pressure [bar]", [s[:Pressure][end] / bar for s in bed]),
        ("air in the gas [-]", [s[:y][i_N2, end] for s in bed]),
        ("temperature [°C]", T .- 273.15),
        ("H₂O loading [mol/kg]", [s[:AdsorbedConcentration][i_H2O, end] / ρ for s in bed]),
        ("relative humidity [-]", [s[:y][i_H2O, end] * s[:Pressure][end] / Mocca.water_vapour_pressure(s[:Temperature][end]) for s in bed]),
    ]
    stages = Mocca.stages(p)
    ends = cumsum([st.duration for st in stages])
    colours = Dict("pressurisation" => "#8a63b8", "adsorption" => "#3d8c47", "blowdown" => "#c48a1f",
        "heating" => "#d9733b", "desorption" => "#3f63c9")
    fig = Figure(size = (900, 640))
    for (k, (label, v)) in enumerate(series)
        ax = Axis(fig[cld(k, 2), mod1(k, 2)], ylabel = label, xlabel = k > 4 ? "time in the cycle [s]" : "")
        for (j, st) in enumerate(stages)
            idx = findall(x -> (j == 1 ? 0.0 : ends[j - 1]) - 1e-6 <= x - result.timesteps[1] / 2 <= ends[j] + 1e-6, t)
            lines!(ax, t[idx], v[idx], color = colours[st.name], linewidth = 2, label = st.name)
        end
        k == 1 && axislegend(ax, position = :rt, framevisible = false, labelsize = 11)
    end
    return fig
end
f_exit = exit_figure(result)

# # 6. Steam or no steam, and the desorption time
#
# Figure 6 of the paper compares cycles with and without steam (its cycles A
# and B) over desorption times from 2000 to 50000 s. A longer desorption
# takes more CO₂ off the bed in each cycle, but makes the cycle longer, so the
# production rate peaks; past the peak the vacuum pump has little left to do,
# while the steam keeps flowing. Each point here is a cyclic steady state.
desorption_times = [3000.0, 6000.0, 10000.0, 15000.0, 20000.0, 30000.0, 50000.0]
designs = [(label = "steam, 25 cm³/s", steam_flow = 25e-6), (label = "no steam", steam_flow = 0.0)]
sweep = map(designs) do d
    map(desorption_times) do t4
        p = Mocca.dac_tvsa_process(; steam_flow = d.steam_flow, durations = (13772.0, 30.0, 704.0, t4))
        dac_performance(Mocca.cyclic_steady_state(p))
    end
end
for (d, rows) in zip(designs, sweep)
    k = argmax([r.production for r in rows])
    println(rpad(d.label, 16), "peak ", round(rows[k].production * 1e5, sigdigits = 3), "×10⁻⁵ kg/h at t₄ = ",
        Int(desorption_times[k]), " s: ", round(rows[k].electrical, sigdigits = 2), " MJ/kg electrical, ",
        round(rows[k].thermal, sigdigits = 3), " MJ/kg thermal")
end

# The production against the specific energies, as in Figure 6 of the paper.
# Along each curve the desorption time rises from 3000 s, the lowest point at
# the right of the electrical energy, to 50000 s, the lowest point at the left:
function sweep_figure(designs, sweep)
    fig = Figure(size = (900, 360))
    for (j, (key, label)) in enumerate(((:electrical, "electrical energy [MJ/kg]"), (:thermal, "thermal energy [MJ/kg]")))
        ax = Axis(fig[1, j], xlabel = label, ylabel = j == 1 ? "CO₂ production [10⁻⁵ kg/h]" : "")
        for (d, rows, colour) in zip(designs, sweep, ("#3f63c9", "#8d8a94"))
            x = [getfield(r, key) for r in rows]
            y = [r.production * 1e5 for r in rows]
            scatterlines!(ax, x, y, color = colour, label = d.label)
        end
        j == 1 && axislegend(ax, position = :rb, framevisible = false)
    end
    return fig
end
f_sweep = sweep_figure(designs, sweep)

# Both cycles produce most at a desorption time of about 10000 s. Without
# steam the peak is 13% lower, for a tenth of the heat; the paper finds 9% less
# CO₂ for 70% less heat. The heat of the cycle without steam, 9 MJ/kg at its
# peak, is close to the paper's 10 MJ/kg; the difference with steam is the
# steam flow, as above. The steam hardly changes the electrical energy, here
# as in the paper.
#
# The [Optimising direct air capture](dac_optimisation.md) example finds these
# designs with gradients instead: the desorption time that produces the most,
# and the steam flow and desorption time that take the least heat for nearly
# as much CO₂.
