# # A column with a cooled wall, across its radius

# A 1D bed has one temperature across its cross-section, and loses heat to
# the wall through a coefficient that lumps the resistance of the whole bed.
# Where the bed conducts heat poorly compared with how fast the wall takes it,
# the centre of a column runs hotter than the gas at the wall during
# adsorption, and a 1D bed cannot show it. This example meshes the bed across
# its radius as well, in rings from the axis to the wall, for a breakthrough
# into a column with a cooled wall, and shows when the 1D bed is enough.

import Mocca
using Mocca.Units
using CairoMakie

# # The column
# CO₂ and N₂ (15/85) on zeolite 13X, as in the [breakthrough](breakthrough.md)
# example, in the column of Haghpanah et al. (2013), here with a wall cooled
# from outside: a large outer coefficient keeps the wall near the ambient
# 25 °C. A bed meshed across its radius needs two coefficients a 1D bed
# lumps together: the effective radial conductivity of the packing, ``K_r``,
# and the coefficient of the film of gas at the wall, ``h_w``. The values
# here are for illustration.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)
adsorbent = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")
column(K_r) = Mocca.Column(; length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    axial_conductivity = 0.0903, outer_heat_transfer = 500.0, radial_conductivity = K_r, wall_heat_transfer = 30.0)

# How much the temperature varies across the radius is set by the Biot
# number of the cross-section, ``\mathrm{Bi} = h_w r_{in} / K_r``: well below
# 1 the cross-section is nearly at one temperature, well above it the wall
# cools only the gas next to it.
biot(c) = c.wall_heat_transfer * c.radius / c.radial_conductivity
biot(column(0.2))

# # A breakthrough on rings
# The same breakthrough on a 1D bed and on 30 cells along the column, each
# split into 6 rings. The wall cools a layer of the bed next to it that is
# thin at this Biot number, so the rings are narrower towards the wall: their
# edges, from the axis, are `radii`. The 1D bed takes the wall coefficient
# that gives the mean temperature of the bed across its radius,
# ``1/h_{in} = 1/h_w + r_{in}/(4 K_r)`` (see [Model equations](@ref "The wall")).
lumped(c) = 1 / (1 / c.wall_heat_transfer + c.radius / (4 * c.radial_conductivity))
radii = 0.1445 .* [0.0, 0.4, 0.65, 0.8, 0.9, 0.96, 1.0]

function breakthrough(c; radial_cells = 1)
    if radial_cells == 1
        c = Mocca.Column(; length = c.length, radius = c.radius, outer_radius = c.outer_radius, porosity = c.porosity,
            axial_conductivity = c.axial_conductivity, outer_heat_transfer = c.outer_heat_transfer,
            inner_heat_transfer = lumped(c))
    end
    fs = Mocca.Flowsheet(gas, feed)
    Mocca.add_bed!(fs, :Bed, adsorbent, c; ncells = 30, radial_cells = radial_cells,
        radii = radial_cells > 1 ? radii : nothing)
    Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))
    p = Mocca.Process(fs, Mocca.Cycle([Mocca.Adsorb(1200.0; velocity = 0.5)]);
        timestepping = Mocca.TimeStepping(max_dt = 10.0, first_dt = 0.5),
        initial_state = Dict(:Bed => (pressure = 1bar, temperature = 298.15, composition = [1e-6, 1 - 1e-6])))
    return p, Mocca.simulate(p)
end

p1, r1 = breakthrough(column(0.2))
p6, r6 = breakthrough(column(0.2); radial_cells = 6)
bed = p6.flowsheet[:Bed]
nothing #hide

# The gas leaving the bed is the rings at the top mixed. On rings it breaks
# through a little later and leaves cooler: the gas next to the cooled wall
# stays cold, and the sorbent there takes up more CO₂ before the front passes.
outlet(p, r) = [Mocca.port_state(p.flowsheet[:Bed], s[:Bed], :top) for s in r.states]
f_outlet = Figure(size = (800, 350))
ax_y = Axis(f_outlet[1, 1], xlabel = "t [s]", ylabel = "CO₂ at the outlet [-]")
ax_T = Axis(f_outlet[1, 2], xlabel = "t [s]", ylabel = "Temperature at the outlet [K]")
for (label, p, r) in (("1D bed", p1, r1), ("30 × 6 rings", p6, r6))
    t = cumsum(r.timesteps)
    lines!(ax_y, t, [o.composition[1] for o in outlet(p, r)], label = label)
    lines!(ax_T, t, [o.temperature for o in outlet(p, r)], label = label)
end
axislegend(ax_y, position = :lt)
f_outlet

# Inside the bed the two differ more. 400 s in, the heat of adsorption has
# warmed the bed behind the front, and the wall has cooled the gas next to it,
# but not the centre:
k = findlast(<=(400.0), cumsum(r6.timesteps))
f_section = Mocca.plot_bed_section(r6.states[k], bed; unit = :Bed, radial_scale = 2)

# Across the radius, at a quarter and half way up the bed:
f_radial = Mocca.plot_radial_profiles(r6.states, r6.timesteps, bed, [200.0, 400.0, 800.0];
    unit = :Bed, positions = [0.25, 0.5])

# The heat leaves through a layer next to the wall, which the narrow outer
# rings resolve. At its most, the centre runs this much hotter than the outer
# ring [K]:
centre_minus_wall(r) = maximum(maximum(s[:Bed][:Temperature][1:30] .- s[:Bed][:Temperature][151:180]) for s in r.states)
centre_minus_wall(r6)

# # When the 1D bed is enough
# As the radial conductivity rises, the Biot number falls and the
# cross-section approaches one temperature. The hottest mean temperature of
# the cross-section half way up the bed, on rings and for the 1D bed with the
# lumped coefficient, and the largest difference between the centre and the
# outer ring:
peak(p, r) = maximum(Mocca.axial_profile(p.flowsheet[:Bed], s[:Bed][:Temperature])[15] for s in r.states)
for K_r in (0.2, 1.0, 10.0)
    c = column(K_r)
    q1, s1 = breakthrough(c)
    q6, s6 = breakthrough(c; radial_cells = 6)
    println("K_r = ", rpad(K_r, 5), " Bi = ", rpad(round(biot(c), digits = 1), 5),
        " rings ", round(peak(q6, s6), digits = 2), " K, 1D ", round(peak(q1, s1), digits = 2),
        " K, centre − wall ", round(centre_minus_wall(s6), digits = 2), " K")
end

# The mean temperature of the 1D bed is within about 1.5 K of that on rings
# even where the centre and the wall differ by 45 K: the lumped coefficient
# gets the mean right, and the performance with it. The rings are for what
# the mean hides, such as the hottest sorbent at the centre, or the cooler and
# more loaded sorbent next to the wall.

# # A cycle on rings
# A bed meshed across its radius works in any process. The four-step VSA of
# Haghpanah et al. (2013) on 20 × 4 cells, to its cyclic steady state, which
# `cyclic_steady_state` finds first for the same bed with one radial cell:
vsa = Mocca.four_stage_vsa_process(ncells = 20, radial_cells = 4, column = column(0.2))
steady = Mocca.cyclic_steady_state(vsa; progress = true)
Mocca.performance(steady)

# The ports are the end rings, so the devices, steps and performance are as
# for any bed. [`write_results`](@ref) adds every cell of the bed in
# `bed_fields.csv`, and the summary the largest temperature difference across
# its radius over the cycle.
