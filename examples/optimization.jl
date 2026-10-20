# # Optimization of a cyclic process

# We optimise the four-stage vacuum swing adsorption (VSA) cycle of
# Haghpanah et al. (2013) at cyclic steady state, twice over:
# 1. maximise the CO₂ recovery, keeping the purity of the CO₂ product at 90%
#    or above;
# 2. from that design, minimise the energy per tonne of CO₂ captured, keeping
#    the purity at 90% and the recovery at 60% or above.
#
# Each evaluation finds the cyclic steady state, starting from the last one,
# and its gradient with respect to the design variables from adjoints (see the
# [gradients](flowsheet_vsa_gradients.md) example). A gradient-based optimiser
# (BFGS within the bounds) then needs a few tens of evaluations, where the
# genetic algorithm of the paper evaluates thousands. For more details about
# the VSA modelling, see the [process from its parts](process_from_parts.md)
# example.

import Mocca
using Mocca.Units
using CairoMakie

# # The process and its design variables
# Four of the decision variables of section 7.1.1 of the paper: the feed
# velocity, the evacuation pressure `P_L`, the blowdown pressure `P_I` and the
# adsorption time. The bounds of the pressures and the velocity are those of
# Table 5 of the paper. The paper also varies the blowdown and evacuation
# times, which we leave out to keep the example short.
#
# Since a duration is among the variables, a stage of a different duration
# keeps the time steps it has now, stretched in proportion, so that the
# performance changes smoothly with the duration.
process = Mocca.four_stage_vsa_process()
Mocca.add_variable!(process, :v_feed, ("adsorption", :velocity); bounds = (0.1, 2.0), description = "interstitial feed velocity [m/s]")
Mocca.add_variable!(process, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.5bar), description = "evacuation pressure")
Mocca.add_variable!(process, :P_I, ("blowdown", :to); bounds = (0.05bar, 0.5bar), description = "blowdown pressure")
Mocca.add_variable!(process, :t_ads, ("adsorption", :duration); bounds = (10.0, 100.0), description = "adsorption time [s]")

# The blowdown must end above the evacuation pressure, which the optimiser
# takes as a linear constraint on the variables, `P_L − P_I ≤ −0.01 bar`:
blowdown_above_evacuation = [(Dict(:P_L => 1.0, :P_I => -1.0), -0.01bar)]

# The design we start from, and a line for each design the optimiser tries:
show_design(values) = join(["$k = " * (k in (:P_L, :P_I) ? "$(round(v / bar, digits = 3)) bar" : "$(round(v, sigdigits = 3))")
    for (k, v) in pairs(values)], ", ")
show_performance(m) = "purity $(round(100m.purity, digits = 1))%, recovery $(round(100m.recovery, digits = 1))%, " *
    "$(round(m.specific_energy, digits = 1)) kWh/t"
report(values, r) = println("  ", show_design(values), ": ", show_performance(r.metrics))

start = Mocca.performance(process; progress = true)
println("start: ", show_performance(start))

# # 1. Maximise the recovery
# The purity target is a penalty, so the optimum may miss it slightly.
max_recovery = Mocca.optimise!(process; maximise = :recovery, at_least = (purity = 0.9,),
    linear_constraints = blowdown_above_evacuation, max_it = 12, log = report)
println("$(max_recovery.evaluations) evaluations: ", show_design(max_recovery.values))
println(show_performance(max_recovery.metrics))

# The evacuation pressure reaches its lower bound within ten evaluations: a
# deeper vacuum takes more CO₂ off the adsorbent, so the recovery rises as `P_L`
# falls. Haghpanah et al. find the same (their Figure 11), and reach 90% purity
# and recovery together only with `P_L` ≤ 0.02 bar. The blowdown pressure then
# sets the purity: blowing down further removes more of the nitrogen, but loses
# more CO₂ with it.

# # 2. Minimise the energy
# The optimiser starts from the design above, where the process was left.
# Each target is a penalty of `weight` (1000 by default) times the square of
# the shortfall, so the objective must be of order one: with the energy in
# units of 100 kWh/t, a shortfall of 1 percentage point of recovery costs as
# much as 10 kWh/t.
min_energy = Mocca.optimise!(process; minimise = :specific_energy, at_least = (purity = 0.9, recovery = 0.6),
    scale = (specific_energy = 100.0,), linear_constraints = blowdown_above_evacuation, max_it = 15, log = report)
println("$(min_energy.evaluations) evaluations: ", show_design(min_energy.values))
println(show_performance(min_energy.metrics))

# # Comparison of the designs
# The process is left at the last optimum. The two optimisations side by side:
designs = ["start" => first(max_recovery.history), "max recovery" => max_recovery, "min energy" => min_energy]
rows = [
    "feed velocity [m/s]" => d -> d.values.v_feed,
    "evacuation P_L [bar]" => d -> d.values.P_L / bar,
    "blowdown P_I [bar]" => d -> d.values.P_I / bar,
    "adsorption time [s]" => d -> d.values.t_ads,
    "purity [%]" => d -> 100d.metrics.purity,
    "recovery [%]" => d -> 100d.metrics.recovery,
    "energy [kWh/t]" => d -> d.metrics.specific_energy,
    "productivity [mol/m³/s]" => d -> d.metrics.productivity,
]
println(rpad("", 24), join(lpad.(first.(designs), 14)))
for (label, f) in rows
    println(rpad(label, 24), join([lpad(round(f(d), sigdigits = 3), 14) for d in last.(designs)]))
end

# Every design tried, on the purity and the energy against the recovery, with
# the targets dashed and each optimum starred:
fig = Figure(size = (1000, 420))
ax_purity = Axis(fig[1, 1], xlabel = "CO₂ recovery [%]", ylabel = "CO₂ purity [%]")
ax_energy = Axis(fig[1, 2], xlabel = "CO₂ recovery [%]", ylabel = "Specific energy [kWh/t]")
hlines!(ax_purity, [90.0], color = :gray, linestyle = :dash)
vlines!(ax_energy, [60.0], color = :gray, linestyle = :dash)
for (r, color, label) in ((max_recovery, :steelblue, "1. max recovery"), (min_energy, :darkorange, "2. min energy"))
    recovery = [100h.metrics.recovery for h in r.history]
    for (ax, y) in ((ax_purity, [100h.metrics.purity for h in r.history]), (ax_energy, [h.metrics.specific_energy for h in r.history]))
        scatterlines!(ax, recovery, y, color = color, label = label, markersize = 6, linewidth = 0.5)
    end
    scatter!(ax_purity, [100r.metrics.recovery], [100r.metrics.purity], color = color, marker = :star5, markersize = 20)
    scatter!(ax_energy, [100r.metrics.recovery], [r.metrics.specific_energy], color = color, marker = :star5, markersize = 20)
end
axislegend(ax_purity, position = :rb)
fig

# Minimising the energy gives up five points of recovery for about a quarter
# less energy per tonne. It does not raise the evacuation pressure, which stays
# at its bound: it blows down less deeply and feeds for longer and faster, so
# that the pump work is spread over more CO₂ captured per cycle. That also
# raises the productivity, so the same plant would need less adsorbent.
