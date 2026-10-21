# Beds and processes of Haghpanah et al. (2013) for the tests. They are built
# from the constants of Mocca 0.1.0 through convert_constants, so they keep the
# axial dispersion at the constants' velocity.

using Mocca
using Mocca.Units

haghpanah_constants(; kwarg...) = Mocca.HaghpanahConstants{Float64}(; kwarg...)

"Model of a bed of the constants `c`."
function haghpanah_bed(c = haghpanah_constants(); ncells = 20, thermal = Mocca.WithWall())
    k = Mocca.convert_constants(c)
    return Mocca.fixed_bed_model(k.gas, k.feed, k.adsorbent, k.column;
        ncells = ncells, thermal = thermal, dispersion_velocity = k.dispersion_velocity)
end

"""
The bed of `c` split along its length into `channels` equal channels side by
side, on `ncells × channels` cells, so that each port is a layer of `channels`
cells. With `porosity`, one value per channel, the channels hold different
amounts of gas and sorbent.
"""
function haghpanah_channel_bed(c = haghpanah_constants(); ncells = 20, channels = 2, thermal = Mocca.Adiabatic(),
        porosity = fill(c.Φ, channels))
    bed = haghpanah_bed(c; ncells = ncells, thermal = thermal)
    d = bed.data_domain
    side = sqrt(π * c.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, channels, 1), (sum(d[:dx]), side, side))
    keys_ = (:r_in, :r_out, :diffusion_coefficient, :thermal_conductivity, :adsorbent_density,
        :adsorbent_heat_capacity, :wall_density, :wall_heat_capacity, :wall_conductivity, :fluid_viscosity,
        :inner_htc, :outer_htc, :ambient_temperature)
    domain = Mocca.mocca_domain(mesh; porosity = repeat(porosity, inner = ncells), permeability = d[:permeability][1],
        (k => first(d[k, Mocca.Unit()]) for k in keys_)...)
    return Mocca.setup_process_model(bed.system, domain)
end

"""
Flowsheet of a bed of `c` with a feed and a product device, and an extract
device if `extract`. With `bed_model`, the bed is that model.
"""
function haghpanah_flowsheet(c; ncells, thermal = Mocca.WithWall(), extract = true, bed_model = nothing)
    k = Mocca.convert_constants(c)
    fs = Flowsheet(k.gas, k.feed)
    add_bed!(fs, :Bed, k.adsorbent, k.column; ncells = ncells, thermal = thermal,
        dispersion_velocity = k.dispersion_velocity)
    isnothing(bed_model) || (Mocca.network(fs).units[:Bed] = bed_model)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    extract && add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    return fs, k
end

"Breakthrough of the feed of `c` into a bed full of nitrogen at the feed pressure, for `duration` seconds."
function haghpanah_breakthrough(c = haghpanah_constants(); duration, max_dt = duration, kwarg...)
    fs, k = haghpanah_flowsheet(c; extract = false, kwarg...)
    return Process(fs, Cycle([Adsorb(duration; velocity = k.velocity)]);
        timestepping = TimeStepping(max_dt = max_dt, first_dt = max_dt),
        initial_state = Dict(:Bed => (pressure = k.p_high, temperature = c.T_feed, wall_temperature = c.T_a)))
end

"The four-stage VSA cycle of `c`, starting from a bed full of nitrogen at `initial_pressure`."
function haghpanah_vsa(c = haghpanah_constants(); durations = (15.0, 15.0, 30.0, 40.0), max_dt = 1.0, first_dt = max_dt,
        initial_pressure = c.p_low, kwarg...)
    fs, k = haghpanah_flowsheet(c; kwarg...)
    t_press, t_ads, t_blow, t_evac = durations
    cycle = Cycle([Pressurise(t_press), Adsorb(t_ads; velocity = k.velocity), Blowdown(t_blow; to = k.p_intermediate),
        Evacuate(t_evac; to = k.p_low)])
    return Process(fs, cycle; ramp_rate = k.λ, timestepping = TimeStepping(max_dt = max_dt, first_dt = first_dt),
        initial_state = Dict(:Bed => (pressure = initial_pressure, temperature = 298.15, wall_temperature = c.T_a)))
end

"States of `unit` in the states of a flowsheet simulation."
unit_states(states, unit = :Bed) = [s[unit] for s in states]
