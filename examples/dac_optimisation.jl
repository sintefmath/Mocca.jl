#=
# Optimising direct air capture

The [Direct air capture by steam-assisted TVSA](dac_stampi_bombelli_2020.md)
example runs the S-TVSA cycle of Stampi-Bombelli et al. (2020) and sweeps its
desorption time with and without steam, as in Figure 6 of the paper. Here we
optimise it instead, with gradients:

1. the desorption time that produces the most CO₂ per hour;
2. the desorption time and steam flow that take the least heat per kg of CO₂,
   while producing at least 95% of that.

The paper's measures are not among Mocca's built-in ones, so we add them as
measures of integrals over the cycle (see [`add_integral!`](@ref
Mocca.add_integral!) and [`add_measure!`](@ref Mocca.add_measure!)): the
production rate in kg/h, and the heat of the jacket and of the steam per kg
of CO₂. The heat of the jacket is a built-in integral, `heating`, which reads
the jacket temperature of each step; the steam is the water that enters
through the steam inlet, raised from water at 20 °C. Each evaluation finds the
cyclic steady state, starting from the last one, and its gradients by
adjoints at the steady state (see [`performance_gradient`](@ref
Mocca.performance_gradient)).

The whole run takes about 13 minutes, so the docs show its results without
running it.
=#

import Mocca
using Mocca.Units

process = Mocca.read_process("stampi_bombelli_2020_stvsa")
gas = process.flowsheet.gas
M = gas.molar_masses
i_H2O = findfirst(==("H2O"), gas.names)

# # The measures
#
# The heat to raise steam at 50 mbar and 95 °C from water at 20 °C, per kg:
# heating the water to its boiling point at 50 mbar, boiling it, and heating
# the steam, as in eq. (19) of the paper.
T_b = 1730.63 / (8.07131 - log10(50mbar / 133.322368)) - 233.426 + 273.15
h_steam = 4186.0 * (T_b - 293.15) + (2.501e6 - 2370.0 * (T_b - 273.15)) + 1864.0 * (368.15 - T_b)

# The moles of water through the steam inlet over a cycle, as an integral, and
# the two measures from the integrals: `s.key_extract` is the CO₂ produced over
# the cycle [mol], `s.duration` the cycle time [s] and `s.heating` the heat of
# the jacket [J].
Mocca.add_integral!(process, :steam, Mocca.StreamObjective(:V_steam; component = i_H2O))
Mocca.add_measure!(process, :production, s -> s.key_extract * M[process.options.key_component] / s.duration * 3600)
Mocca.add_measure!(process, :thermal_energy, s -> (s.heating + s.steam * M[i_H2O] * h_steam) / (s.key_extract * M[process.options.key_component]) / 1e6)

# At the base case of the paper, with 30000 s of desorption, the cycle produces
# 2.19×10⁻⁵ kg of CO₂ per hour, for 218 MJ of heat per kg, nearly all for the
# steam:
base = Mocca.performance(process)
println("production ", round(base.production * 1e5, sigdigits = 3), "×10⁻⁵ kg/h, heat ", round(base.thermal_energy, sigdigits = 3), " MJ/kg")

# # 1. The most CO₂ per hour
#
# A longer desorption takes more CO₂ off the bed in each cycle, but makes the
# cycle longer. The gradient of the production rate with respect to the
# desorption time finds the balance:
Mocca.add_variable!(process, :t_des, ("desorption", :duration); bounds = (3000.0, 30000.0),
    description = "desorption time [s]")
show_design(v) = join(["$k = " * (k == :steam_flow ? "$(round(x * 1e6, sigdigits = 3)) cm³/s" : "$(round(Int, x)) s") for (k, x) in pairs(v)], ", ")
show_metrics(m) = "$(round(m.production * 1e5, sigdigits = 4))×10⁻⁵ kg/h, $(round(m.thermal_energy, sigdigits = 4)) MJ/kg"
report(v, r) = println("  ", show_design(v), ": ", show_metrics(r.metrics))
most = Mocca.optimise!(process; maximise = :production, scale = (production = 1e-5,), max_it = 10, log = report)
println("$(most.evaluations) evaluations: ", show_design(most.values), ", ", show_metrics(most.metrics))

# It takes 9 evaluations to settle at a desorption time of 10471 s, which
# produces 3.22×10⁻⁵ kg/h, 47% more than the base case, for 97 MJ/kg. The
# sweep of the other example has its peak there too.

# # 2. The least heat
#
# Nearly all the heat goes to the steam, and the steam helps desorption less
# than its cost: the paper finds that without it the production peaks 9% lower
# for 70% less heat. Here the steam flow is a variable too, between 2 and
# 25 cm³/s, and the heat per kg of CO₂ is minimised, keeping the production at
# 95% of its most or above:
Mocca.add_variable!(process, :steam_flow, ("desorption", :flow_rate); bounds = (2e-6, 25e-6),
    description = "steam flow at 50 mbar and 95 °C [m³/s]")
least = Mocca.optimise!(process; minimise = :thermal_energy, at_least = (production = 0.95 * most.metrics.production,),
    scale = (thermal_energy = 10.0, production = 1e-5), max_it = 12, log = report)
println("$(least.evaluations) evaluations, feasible: $(least.feasible): ", show_design(least.values), ", ", show_metrics(least.metrics))

# After 34 evaluations the steam is down to 4.0 cm³/s, with a slightly longer
# desorption of 11295 s to make up for it, and the heat to 24 MJ/kg: 76% less
# than at the most productive design, for the 5% less CO₂ allowed. The
# production target holds the steam above its lower bound: at 2 cm³/s the
# production falls to 3.0×10⁻⁵ kg/h, just under it.

# # The designs
designs = ["base case" => (values = (t_des = 30000.0, steam_flow = 25e-6), metrics = base),
    "most CO₂" => (values = (t_des = most.values.t_des, steam_flow = 25e-6), metrics = most.metrics),
    "least heat" => least]
println(rpad("", 26), join(lpad.(first.(designs), 14)))
for (label, f) in ["desorption time [s]" => d -> round(Int, d.values.t_des),
        "steam [cm³/s]" => d -> round(d.values.steam_flow * 1e6, sigdigits = 3),
        "production [1e-5 kg/h]" => d -> round(d.metrics.production * 1e5, sigdigits = 3),
        "heat [MJ/kg]" => d -> round(d.metrics.thermal_energy, sigdigits = 3),
        "recovery [%]" => d -> round(100 * d.metrics.recovery, sigdigits = 3)]
    println(rpad(label, 26), join([lpad(f(last(d)), 14) for d in designs]))
end
