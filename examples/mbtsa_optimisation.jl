#=
# Optimising the moving bed TSA

The [Moving bed TSA](mbtsa_mondino_2022.md) example builds the MBTSA process
of Mondino et al. (2022) and finds how its capture rate and its work per tonne
of CO₂ respond to its settings. Here we optimise it: the least work per tonne
of CO₂ that keeps the purity at 95% and the capture rate at 90%, the targets of
the paper, and then what a more productive plant costs in work.

The work counts the electrical work of the compressor, the vacuum pump and the
blowers, and the work equivalent to the heat, each bit of heat weighted by its
Carnot factor against 25 °C (`specific_equivalent_work`; see
[`performance_integrals`](@ref Mocca.performance_integrals)). Each evaluation
runs the plant to steady state, starting from the last one, and takes the
gradients by adjoints at the steady state (see
[`performance_gradient`](@ref Mocca.performance_gradient)).

The process comes from Mocca's library, on a grid of 50 cells per metre,
coarser than the 200 of that example, to keep the run short: its capture rate
is 1.4 points lower (91.7% against 93.1%), which matters little when designs
are compared with one another.
=#

import Mocca
using Mocca.Units

process = Mocca.read_process("mondino_2022_mbtsa")
gas = process.flowsheet.gas
feed = process.flowsheet.feed

# # 1. The design variables
# Four settings, within bounds that keep the plant in the regime of the paper:
#
# * the purge, as a percentage of the flue gas by mass (55.9 kg/s): the purge
#   rate is a setting of the device `V_purge`, in mol/s, and a
#   [`ScaledPlace`](@ref Mocca.ScaledPlace) sets it in proportion. Below about
#   14% no purge gas is left to strip CO₂ from the hot sorbent and the capture
#   rate falls by four points (see the moving bed example), so the bound is
#   14.5%;
# * the temperature of the hot water heating the desorption section;
# * the sorbent flow, one variable for the flow of all five sections;
# * the height of the desorption section, a
#   [`GeometryPlace`](@ref Mocca.GeometryPlace): the section is rebuilt on as
#   many cells.
#
# The sorbent-side heat transfer coefficients stay at the values computed from
# the sorbent velocity of the paper's design.
M_purge = sum(feed.light_product .* gas.molar_masses)  # kg/mol
mol_per_percent = 55.9 / 100 / M_purge
Mocca.add_variable!(process, :purge, Mocca.ScaledPlace(Mocca.SettingPlace(:V_purge, :rate), mol_per_percent);
    bounds = (14.5, 20.0), description = "purge [% of the flue gas by mass]")
Mocca.add_variable!(process, :hot_water, (:parameter, :Desorption, :HXFluidInletTemperature);
    bounds = (celsius(167), celsius(197)), description = "hot water temperature [K]")
sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
Mocca.add_variable!(process, :sorbent_flow, ((:parameter, s, :SorbentMassFlow) for s in sections)...;
    bounds = (550.0, 750.0), description = "sorbent flow [kg/s]")
Mocca.add_variable!(process, :desorption_height, (:geometry, :Desorption, :length);
    bounds = (0.6, 1.0), description = "height of the desorption section [m]")

show_design(v) = "purge $(round(v.purge, digits = 2))%, hot water $(round(v.hot_water - 273.15, digits = 1)) °C, " *
    "sorbent $(round(v.sorbent_flow, digits = 1)) kg/s, desorption $(round(v.desorption_height, digits = 3)) m"
show_performance(m) = "purity $(round(100m.purity, digits = 2))%, capture $(round(100m.recovery, digits = 2))%, " *
    "work $(round(3.6e-3 * m.specific_equivalent_work, digits = 3)) MJ/kg, productivity $(round(m.productivity, sigdigits = 4)) mol/(m³·s)"
design(p) = NamedTuple{Tuple(Mocca.variable_names(p))}(Tuple(Mocca.variable_values(p)))
start = Mocca.performance(process)
println("start: ", show_design(design(process)))
println("       ", show_performance(start))

# # 2. A look around
# Before following gradients, a few designs spread over the box of the
# bounds, each run from the initial state (see [`explore`](@ref Mocca.explore)),
# show where the targets can be met:
points = Mocca.explore(process; samples = 6, log = (v, m) -> println("  ", show_design(v), ": ",
    isnothing(m) ? "failed" : show_performance(m)))
println(count(p -> p.settled, points), " of ", length(points), " designs settled")

# The paper's design gives 99.1% purity and 91.7% capture, more than the
# targets ask, for 1.83 MJ/kg. Cooler water saves work, as the moving bed
# example found: the designs with water at 170–177 °C need about 1.6 MJ/kg.
# But less sorbent loses capture, 79–87% at 558–598 kg/s, and cool water
# loses purity: the targets cut the box, and the least work lies on their
# edge.

# # 3. The least work for 95% purity and 90% capture
# The work in units of 1 MJ/kg (278 kWh/t), so that the objective is of order
# one. Each BFGS run takes up to ten iterations, and the targets are met by an
# augmented Lagrangian (see [`optimise!`](@ref Mocca.optimise!)).
report(v, r) = println("  ", show_design(v), ": ", show_performance(r.metrics))
best = Mocca.optimise!(process; minimise = :specific_equivalent_work, at_least = (purity = 0.95, recovery = 0.9),
    scale = (specific_equivalent_work = 277.8,), max_it = 10, log = report)
println(best.evaluations, " evaluations; ", best.feasible ? "the targets are met" : "the targets are not met")
println(show_design(best.values))
println(show_performance(best.metrics))

# The optimum, about 50 evaluations, meets both targets exactly and needs
# 1.53 MJ/kg, a sixth less than the paper's design. It lowers the hot water
# to 168 °C, close to its bound, and makes up for the slower desorption with
# a taller desorption section, 0.94 m, and a little more purge, 16.5%: heat
# given at a lower temperature is worth less work, and the taller section
# gives the sorbent the time to release its CO₂ at that temperature.

# # 4. Work against productivity
# A more productive plant captures more CO₂ per cubic metre of sorbent bed, so
# it is smaller for the same flue gas. The least work for productivities
# above the optimum's, each found from the last ([`pareto_front!`](@ref
# Mocca.pareto_front!)), traces how much work a smaller plant costs:
p0 = best.metrics.productivity
front = Mocca.pareto_front!(process; minimise = :specific_equivalent_work, sweep = :productivity => p0 .* [1.05, 1.1, 1.15],
    at_least = (purity = 0.95, recovery = 0.9), scale = (specific_equivalent_work = 277.8, productivity = p0), max_it = 8)
println(rpad("productivity", 14), rpad("work [MJ/kg]", 14), "design")
println(rpad(round(p0, sigdigits = 4), 14), rpad(round(3.6e-3 * best.metrics.specific_equivalent_work, digits = 3), 14), show_design(best.values))
for pt in front
    println(rpad(round(pt.metrics.productivity, sigdigits = 4), 14), rpad(round(3.6e-3 * pt.metrics.specific_equivalent_work, digits = 3), 14),
        show_design(pt.values), pt.feasible ? "" : "  (targets not met)")
end

# A plant 10% more productive needs a desorption section of 0.69 m instead of
# 0.94 m and water at 171 °C, for about 1.5% more work, and one 15% more
# productive the shortest section allowed, 0.6 m, with water at 174 °C and a
# little more sorbent, for about 5% more: the trade-off is shallow at first,
# so a somewhat smaller plant costs little. These optima hold for this model
# and grid; the sorbent-side heat transfer coefficients were not recomputed
# for the new sorbent flows, and the hot water lies near its lower bound.
