# Ready-made processes from the literature, built on the user interface layer.
# Each is a short composition of a gas, a feed, library adsorbents, columns,
# devices with roles and a cycle of step templates, and serves as a template
# for new processes.

# ------------------------------------------------------------------------------
# Haghpanah et al. (2013): CO2/N2 on zeolite 13X
# ------------------------------------------------------------------------------

"""
    haghpanah_2013_gas()

CO2/N2 with the properties used by Haghpanah et al. (2013), Ind. Eng. Chem.
Res. 52, 4249.
"""
haghpanah_2013_gas() = Gas(["CO2", "N2"]; molar_masses = [44.01e-3, 28e-3], heat_capacities = [697.5687, 1096.4],
    viscosity = 1.72e-5, density = 1.22638310956, diffusivity = 1.6e-5)

"""
    haghpanah_2013_column(; length = 1.0)

The column of Haghpanah et al. (2013).
"""
haghpanah_2013_column(; length = 1.0) = Column(; length = length, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    axial_conductivity = 0.0903, wall = Wall(density = 7800.0, heat_capacity = 502.0, conductivity = 16.0),
    inner_heat_transfer = 8.6, outer_heat_transfer = 2.5, ambient_temperature = 298.15)

"""
    four_stage_vsa_process(; ncells = 30, durations = (15.0, 15.0, 30.0, 40.0),
        P_H = 1e5, P_I = 0.2e5, P_L = 0.1e5, velocity = 1.0, ramp_rate = 0.5,
        timestepping = TimeStepping(), name = "four-stage VSA")

The single-bed, four-stage VSA cycle of Haghpanah et al. (2013): CO2/N2
(15/85) on zeolite 13X, with pressurisation, adsorption at interstitial
`velocity` [m/s], co-current blowdown to `P_I` and counter-current evacuation
to `P_L` [Pa].

    feed ──V_feed──┐                 ┌──V_product── raffinate
                   bottom  Bed  top
    extract ──V_vacuum──┘

As in the paper, vacuum pumps discharge at, and the feed compressor takes gas
from, 1 bar.
"""
function four_stage_vsa_process(; ncells = 30, durations = (15.0, 15.0, 30.0, 40.0), P_H = 1e5, P_I = 0.2e5, P_L = 0.1e5,
        velocity = 1.0, ramp_rate = 0.5, timestepping = TimeStepping(), name = "four-stage VSA", initial_state = Dict{Symbol, Any}(),
        adsorbent = load_adsorbent("zeolite_13x_haghpanah_2013"), column = haghpanah_2013_column(), gas = haghpanah_2013_gas())
    feed = Feed(gas; composition = [0.15, 0.85], temperature = 298.15, pressure = P_H)
    fs = Flowsheet(gas, feed)
    add_bed!(fs, :Bed, adsorbent, column; ncells = ncells)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    t_press, t_ads, t_blow, t_evac = durations
    cycle = Cycle([
        Pressurise(t_press),
        Adsorb(t_ads; velocity = velocity),
        Blowdown(t_blow; to = P_I),
        Evacuate(t_evac; to = P_L),
    ])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state, ambient_pressure = 1e5)
end

"""
    two_bed_vsa_process(; ncells = 30, feed_time = 15.0, pressurisation_time = 15.0,
        equalisation_time = 10.0, equalisation_conductance = 1e-6, P_H = 1e5, P_I = 0.2e5,
        P_L = 0.1e5, velocity = 1.0, ramp_rate = 0.5, timestepping = TimeStepping())

Two beds, A and B, running the four-stage VSA cycle half a cycle apart, with
a pressure equalisation valve joining their tops; see
`two_bed_vsa_flowsheet` (deprecated) for the cycle table. After equalising, a
bed is taken to be midway between the high and low pressures, where its
pressurisation starts; evacuation ramps from `P_I`.
"""
function two_bed_vsa_process(; ncells = 30, feed_time = 15.0, pressurisation_time = 15.0, equalisation_time = 10.0,
        equalisation_conductance = 1e-6, P_H = 1e5, P_I = 0.2e5, P_L = 0.1e5, velocity = 1.0, ramp_rate = 0.5,
        timestepping = TimeStepping(), name = "two-bed VSA", initial_state = Dict{Symbol, Any}(),
        adsorbent = load_adsorbent("zeolite_13x_haghpanah_2013"), column = haghpanah_2013_column(), gas = haghpanah_2013_gas())
    feed = Feed(gas; composition = [0.15, 0.85], temperature = 298.15, pressure = P_H)
    fs = Flowsheet(gas, feed)
    for b in (:A, :B)
        add_bed!(fs, b, adsorbent, column; ncells = ncells)
    end
    for b in (:A, :B)
        add_device!(fs, Symbol(:V_feed_, b), :feed; at = (b, :bottom))
        add_device!(fs, Symbol(:V_product_, b), :product; at = (b, :top))
        add_device!(fs, Symbol(:V_vacuum_, b), :extract; at = (b, :bottom))
    end
    add_device!(fs, :V_eq, :link; from = (:A, :top), to = (:B, :top))
    P_eq = (P_H + P_L) / 2
    tp, tf, te = pressurisation_time, feed_time, equalisation_time
    evacuate(b, t) = Evacuate(t; bed = b, from = P_I, to = P_L)
    cycle = Cycle([
        Concurrent(Pressurise(tp; bed = :A, from = P_eq), evacuate(:B, tp)),
        Concurrent(Adsorb(tf; bed = :A, velocity = velocity), evacuate(:B, tf)),
        Equalise(te; beds = (:A, :B), conductance = equalisation_conductance),
        Concurrent(Pressurise(tp; bed = :B, from = P_eq), evacuate(:A, tp)),
        Concurrent(Adsorb(tf; bed = :B, velocity = velocity), evacuate(:A, tf)),
        Equalise(te; beds = (:B, :A), conductance = equalisation_conductance),
    ])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state, ambient_pressure = 1e5)
end

# ------------------------------------------------------------------------------
# Krishnamurthy et al. (2014): wet flue gas
# ------------------------------------------------------------------------------

"""
    krishnamurthy_2014_gas()

Wet flue gas, CO2/N2/H2O, with the properties of Table S2 of Krishnamurthy
et al. (2014), Ind. Eng. Chem. Res. 53, 14462: one gas heat capacity for the
mixture, and a diffusivity per component (CO2–N2 for CO2 and N2, N2–H2O for
water).
"""
krishnamurthy_2014_gas() = Gas(["CO2", "N2", "H2O"]; molar_masses = [44.01e-3, 28e-3, 18.015e-3], heat_capacities = fill(1030.0, 3),
    viscosity = 1.7e-5, density = 1.22638310956, diffusivity = [1.6e-5, 1.6e-5, 2.6e-5])

"""
    krishnamurthy_2014_column(; length = 1.0)

The column of Table S2 of Krishnamurthy et al. (2014).
"""
krishnamurthy_2014_column(; length = 1.0) = Column(; length = length, radius = 0.15, outer_radius = 0.16, porosity = 0.37,
    axial_conductivity = 0.09, wall = Wall(density = 7800.0, heat_capacity = 500.0, conductivity = 16.0),
    inner_heat_transfer = 8.6, outer_heat_transfer = 2.5, ambient_temperature = 298.15)

_wet_feed(gas, P) = Feed(gas; composition = WET_FLUE_GAS_FEED, temperature = 298.15, pressure = P, light_product = (N2 = 1.0,))

"""
    lpp_vsa_process(; ncells = 30, velocity = 0.45, t_lpp = 20.0, t_ads = 78.4, t_bd = 37.9,
        t_evac = 111.6, P_H = 101325.0, P_I = 0.08e5, P_L = 0.03e5, ramp_rate = 0.5,
        timestepping = TimeStepping())

The single-bed VSA cycle with light product pressurisation (LPP) on zeolite
13X for wet flue gas, of section 3 of Krishnamurthy et al. (2014), at the
operating conditions of its section 3.1. The light product is nitrogen.

    feed ──V_feed──┐                    ┌──V_product── raffinate
                   bottom   Bed   top ──┤
    CO2 ◄─V_vacuum─┘                    └──V_lpp────── light product
"""
function lpp_vsa_process(; ncells = 30, velocity = 0.45, t_lpp = 20.0, t_ads = 78.4, t_bd = 37.9, t_evac = 111.6,
        P_H = 101325.0, P_I = 0.08e5, P_L = 0.03e5, ramp_rate = 0.5, timestepping = TimeStepping(), name = "LPP VSA, wet flue gas", initial_state = Dict{Symbol, Any}(),
        adsorbent = load_adsorbent("zeolite_13x_krishnamurthy_2014"), column = krishnamurthy_2014_column(), gas = krishnamurthy_2014_gas())
    fs = Flowsheet(gas, _wet_feed(gas, P_H))
    add_bed!(fs, :Bed, adsorbent, column; ncells = ncells)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_lpp, :light_inlet; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    cycle = Cycle([
        Pressurise(t_lpp; source = :light_product),
        Adsorb(t_ads; velocity = velocity),
        Blowdown(t_bd; to = P_I),
        Evacuate(t_evac; to = P_L),
    ])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state)
end

"""
    dual_adsorbent_vsa_process(; ncells = 30, silica_gel_length = 0.41, zeolite_length = 1.0,
        velocity = 0.70, t_press = 20.0, t_ads = 46.20,
        silica_gel = (t_bd = 42.21, t_evac = 46.02, P_I = 0.48e5, P_L = 0.30e5),
        zeolite = (t_bd = 56.30, t_evac = 101.20, P_I = 0.07e5, P_L = 0.03e5),
        P_H = 101325.0, ramp_rate = 0.5, timestepping = TimeStepping())

The dual-adsorbent, two-bed VSA process of section 4 of Krishnamurthy et al.
(2014): a silica gel bed that takes up the water, in series with a zeolite
13X bed that concentrates the CO2, at the minimum-energy conditions of its
section 4.3. The beds' blowdown and evacuation times differ, so the cycle
runs one sequence per bed side by side.

    feed ──V_feed──┐                      ┌──V_link──┐                   ┌──V_product── raffinate
                   bottom  SilicaGel  top ┘          bottom  Zeolite  top┤
    H2O ◄─V_waste──┘                                 │                   └──V_lpp────── light product
                                          CO2 ◄─V_vacuum┘
"""
function dual_adsorbent_vsa_process(; ncells = 30, silica_gel_length = 0.41, zeolite_length = 1.0, velocity = 0.70,
        t_press = 20.0, t_ads = 46.20,
        silica_gel = (t_bd = 42.21, t_evac = 46.02, P_I = 0.48e5, P_L = 0.30e5),
        zeolite = (t_bd = 56.30, t_evac = 101.20, P_I = 0.07e5, P_L = 0.03e5),
        P_H = 101325.0, ramp_rate = 0.5, timestepping = TimeStepping(), name = "dual-adsorbent VSA, wet flue gas", initial_state = Dict{Symbol, Any}(),
        gas = krishnamurthy_2014_gas())
    fs = Flowsheet(gas, _wet_feed(gas, P_H))
    add_bed!(fs, :SilicaGel, load_adsorbent("silica_gel_krishnamurthy_2014"), krishnamurthy_2014_column(; length = silica_gel_length); ncells = ncells)
    add_bed!(fs, :Zeolite, load_adsorbent("zeolite_13x_krishnamurthy_2014"), krishnamurthy_2014_column(; length = zeolite_length); ncells = ncells)
    add_device!(fs, :V_feed, :feed; at = (:SilicaGel, :bottom))
    add_device!(fs, :V_waste, :waste; at = (:SilicaGel, :bottom))
    add_device!(fs, :V_link, :link; from = (:SilicaGel, :top), to = (:Zeolite, :bottom))
    add_device!(fs, :V_product, :product; at = (:Zeolite, :top))
    add_device!(fs, :V_lpp, :light_inlet; at = (:Zeolite, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Zeolite, :bottom))
    sg, x = silica_gel, zeolite
    cycle = Cycle(
        "silica gel" => [
            Pressurise(t_press; bed = :SilicaGel),
            Adsorb(t_ads; bed = :SilicaGel, velocity = velocity),
            Blowdown(sg.t_bd; bed = :SilicaGel, to = sg.P_I, through = :bottom),
            Evacuate(sg.t_evac; bed = :SilicaGel, to = sg.P_L)],
        "13X" => [
            Pressurise(t_press; bed = :Zeolite, source = :light_product),
            Adsorb(t_ads; bed = :Zeolite),
            Blowdown(x.t_bd; bed = :Zeolite, to = x.P_I),
            Evacuate(x.t_evac; bed = :Zeolite, to = x.P_L)])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state)
end

# ------------------------------------------------------------------------------
# From the old constants structs
# ------------------------------------------------------------------------------

"""
    convert_constants(constants) → NamedTuple

The [`Gas`](@ref), [`Feed`](@ref), [`Adsorbent`](@ref) and [`Column`](@ref)
described by a `HaghpanahConstants` or `adsorptionConstants` struct, and its
operating conditions: `p_high`, `p_intermediate` and `p_low` [Pa], the
superficial feed velocity `v_feed` and the interstitial `velocity` [m/s], the
ramp rate `λ` [1/s] and the `dispersion_velocity` [m/s] the struct used. The
adsorbent counts the heat of sorption of its first component only, as Mocca
0.1.0 did.
"""
function convert_constants(c::ConstantsStruct)
    names = c isa HaghpanahConstants ? ["CO2", "N2"] : collect(String, c.component_names)
    M = c isa HaghpanahConstants ? [c.molecularMassOfCO2, c.molecularMassOfN2] : collect(c.molecular_masses)
    gas = Gas(names; molar_masses = M, heat_capacities = collect(c.C_pg), viscosity = c.fluid_viscosity, density = c.ρ_g, diffusivity = c.D_m)
    iso = dual_site_langmuir((n => (qsb = c.qsbi[i], b0 = c.b0[i], ΔUb = c.ΔUbi[i], qsd = c.qsdi[i], d0 = c.d0[i], ΔUd = c.ΔUdi[i]) for (i, n) in enumerate(names))...;
        reference_temperature = c.T0)
    adsorbent = Adsorbent(""; isotherm = iso, density = c.ρ_s, heat_capacity = c.C_ps, particle_diameter = c.d_p,
        particle_porosity = c.ϵ_p, tortuosity = c.τ, adsorbed_heat_capacity = collect(c.C_pa), sorption_heat = :first_component)
    column = Column(; length = c.L, radius = c.r_in, outer_radius = c.r_out, porosity = c.Φ, axial_conductivity = c.K_z,
        wall = Wall(density = c.ρ_w, heat_capacity = c.C_pw, conductivity = c.K_w),
        inner_heat_transfer = c.h_in, outer_heat_transfer = c.h_out, ambient_temperature = c.T_a)
    feed = Feed(gas; composition = collect(c.y_feed), temperature = c.T_feed, pressure = c.p_high)
    return (gas = gas, feed = feed, adsorbent = adsorbent, column = column, p_high = c.p_high, p_intermediate = c.p_intermediate,
        p_low = c.p_low, v_feed = c.v_feed, velocity = c.v_feed / c.Φ, λ = c.λ, dispersion_velocity = c.V0_inter)
end
