# Ready-made processes from the literature, built on the user interface layer.
# Each is a short composition of a gas, a feed, library adsorbents, columns,
# devices with roles and a cycle of step templates, and serves as a template
# for new processes.

# ------------------------------------------------------------------------------
# Haghpanah et al. (2013): CO2/N2 on zeolite 13X
# ------------------------------------------------------------------------------

"""
    haghpanah_2013_gas(; kwarg...)

CO2/N2 with the properties used by Haghpanah et al. (2013), Ind. Eng. Chem.
Res. 52, 4249. Keywords go to [`Gas`](@ref), such as `transport =
:chapman_enskog`.
"""
haghpanah_2013_gas(; kwarg...) = Gas(["CO2", "N2"]; molar_masses = [44.01e-3, 28e-3], heat_capacities = [697.5687, 1096.4],
    viscosity = 1.72e-5, density = 1.22638310956, diffusivity = 1.6e-5, kwarg...)

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
from, 1 bar, with a heat capacity ratio of 1.4.
"""
function four_stage_vsa_process(; ncells = 30, durations = (15.0, 15.0, 30.0, 40.0), P_H = 1e5, P_I = 0.2e5, P_L = 0.1e5,
        velocity = 1.0, ramp_rate = 0.5, timestepping = TimeStepping(), name = "four-stage VSA", initial_state = Dict{Symbol, Any}(),
        adsorbent = load_adsorbent("zeolite_13x_haghpanah_2013"), column = haghpanah_2013_column(), gas = haghpanah_2013_gas())
    feed = Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = P_H)
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
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state, ambient_pressure = 1e5,
        heat_capacity_ratio = 1.4)
end

"""
    two_bed_vsa_process(; ncells = 30, feed_time = 15.0, pressurisation_time = 15.0,
        equalisation_time = 10.0, equalisation_conductance = 1e-6, P_H = 1e5, P_L = 0.1e5,
        velocity = 1.0, ramp_rate = 0.5, timestepping = TimeStepping())

Two beds, A and B, running the four-stage VSA cycle half a cycle apart, with
a pressure equalisation valve joining their tops:

    feed ──V_feed_A──┐              ┌──V_product_A── raffinate
                     bottom  A  top ┤
    extract ─V_vacuum_A┘            V_eq
    feed ──V_feed_B──┐              │
                     bottom  B  top ┤
    extract ─V_vacuum_B┘            └──V_product_B── raffinate

After a bed finishes adsorption it equalises with the other, evacuated bed,
which saves feed gas and compression work:

| Stage | Bed A                  | Bed B                  |
|-------|------------------------|------------------------|
| 1     | pressurisation (feed)  | evacuation             |
| 2     | adsorption             | evacuation             |
| 3     | equalisation (gives)   | equalisation (takes)   |
| 4     | evacuation             | pressurisation (feed)  |
| 5     | evacuation             | adsorption             |
| 6     | equalisation (takes)   | equalisation (gives)   |

After equalising, a bed is taken to be midway between the high and low
pressures, where its pressurisation or evacuation starts. A bed evacuates
through two stages, along one ramp to `P_L`.
"""
function two_bed_vsa_process(; ncells = 30, feed_time = 15.0, pressurisation_time = 15.0, equalisation_time = 10.0,
        equalisation_conductance = 1e-6, P_H = 1e5, P_L = 0.1e5, velocity = 1.0, ramp_rate = 0.5,
        timestepping = TimeStepping(), name = "two-bed VSA", initial_state = Dict{Symbol, Any}(),
        adsorbent = load_adsorbent("zeolite_13x_haghpanah_2013"), column = haghpanah_2013_column(), gas = haghpanah_2013_gas())
    feed = Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = P_H)
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
    # A bed evacuates through two stages. The second carries on along the ramp
    # of the first, from where it left the bed: the pressure of the ramp from
    # P_eq after the first stage
    P_mid = evaluate_profile(ExponentialRamp(P_eq, P_L, ramp_rate), tp)
    evacuate(b, t, from) = Evacuate(t; bed = b, from = from, to = P_L)
    cycle = Cycle([
        Concurrent(Pressurise(tp; bed = :A, from = P_eq), evacuate(:B, tp, P_eq)),
        Concurrent(Adsorb(tf; bed = :A, velocity = velocity), evacuate(:B, tf, P_mid)),
        Equalise(te; beds = (:A, :B), conductance = equalisation_conductance),
        Concurrent(Pressurise(tp; bed = :B, from = P_eq), evacuate(:A, tp, P_eq)),
        Concurrent(Adsorb(tf; bed = :B, velocity = velocity), evacuate(:A, tf, P_mid)),
        Equalise(te; beds = (:B, :A), conductance = equalisation_conductance),
    ])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state, ambient_pressure = 1e5,
        heat_capacity_ratio = 1.4)
end

# ------------------------------------------------------------------------------
# Krishnamurthy et al. (2014): wet flue gas
# ------------------------------------------------------------------------------

"""
    krishnamurthy_2014_gas(; kwarg...)

Wet flue gas, CO2/N2/H2O, with the properties of Table S2 of Krishnamurthy
et al. (2014), Ind. Eng. Chem. Res. 53, 14462: one gas heat capacity for the
mixture, and a diffusivity per component (CO2–N2 for CO2 and N2, N2–H2O for
water). Keywords go to [`Gas`](@ref).
"""
krishnamurthy_2014_gas(; kwarg...) = Gas(["CO2", "N2", "H2O"]; molar_masses = [44.01e-3, 28e-3, 18.015e-3], heat_capacities = fill(1030.0, 3),
    viscosity = 1.7e-5, density = 1.22638310956, diffusivity = [1.6e-5, 1.6e-5, 2.6e-5], kwarg...)

"""
    krishnamurthy_2014_column(; length = 1.0)

The column of Table S2 of Krishnamurthy et al. (2014).
"""
krishnamurthy_2014_column(; length = 1.0) = Column(; length = length, radius = 0.15, outer_radius = 0.16, porosity = 0.37,
    axial_conductivity = 0.09, wall = Wall(density = 7800.0, heat_capacity = 500.0, conductivity = 16.0),
    inner_heat_transfer = 8.6, outer_heat_transfer = 2.5, ambient_temperature = 298.15)

"Feed composition of the wet flue gas: 15% CO2, 82% N2 and 3% H2O."
const WET_FLUE_GAS_FEED = [0.15, 0.82, 0.03]

# By name, so that the gas may list its components in any order
_wet_feed(gas, P) = Feed(gas; composition = Dict(zip(("CO2", "N2", "H2O"), WET_FLUE_GAS_FEED)), temperature = 298.15, pressure = P,
    light_product = (N2 = 1.0,))

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
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state,
        heat_capacity_ratio = 1.4)
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
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state,
        heat_capacity_ratio = 1.4)
end

# ------------------------------------------------------------------------------
# Stampi-Bombelli et al. (2020): direct air capture
# ------------------------------------------------------------------------------

"""
    stampi_bombelli_2020_gas()

Humid air as CO2, H2O and N2, as in Stampi-Bombelli et al. (2020),
Adsorption 26, 1183, which models the oxygen and argon of air as nitrogen.
"""
stampi_bombelli_2020_gas() = Gas(["CO2", "H2O", "N2"])

"""
    stampi_bombelli_2020_column()

The column of Table 1 of Stampi-Bombelli et al. (2020): a bed 1 cm long and
8 cm across of APDES-NFC at 55.4 kg/m³, in a wall with a heat capacity of
4 MJ/(m³·K) and a heating jacket outside it. The porosity is the paper's
total void fraction, `1 − 55.4/1590`, with the density of the material in
[`load_adsorbent`](@ref)`("apdes_nfc_stampi_bombelli_2020")`. Axial
conduction in the bed and in the wall is neglected, as in the paper.
"""
stampi_bombelli_2020_column() = Column(; length = 0.01, radius = 0.04, outer_radius = 0.041, porosity = 1 - 55.4 / 1590,
    axial_conductivity = 0.0, wall = Wall(density = 8000.0, heat_capacity = 500.0, conductivity = 0.0),
    inner_heat_transfer = 3.0, outer_heat_transfer = 26.0, ambient_temperature = 293.15)

"""
    dac_tvsa_process(; ncells = 10, air_flow = 50e-6, steam_flow = 25e-6, P_L = 0.05e5,
        T_H = 368.15, durations = (13772.0, 30.0, 704.0, 30000.0), t_press = 30.0,
        ramp_rate = 0.5, timestepping = TimeStepping(max_dt = 100.0, first_dt = 0.01),
        permeability = nothing, dispersion = 0.0)

Direct air capture by steam-assisted temperature–vacuum swing adsorption
(S-TVSA) on an amine sorbent, after Stampi-Bombelli et al. (2020), Adsorption
26, 1183, at its base case (Table 5). Humid air at 20 °C, with 400 ppm CO2 and
50% relative humidity, flows at `air_flow` [m³/s] through a short bed of
APDES-NFC. The cycle has the paper's four steps, after a short
pressurisation:

0. pressurisation with air through the feed end, over `t_press` [s]. The
   paper instead opens the bed to 1 bar at the start of adsorption, which is
   too abrupt a jump for the solver;
1. adsorption at 1 bar, with the jacket at ambient temperature;
2. blowdown through the vent to `P_L` [Pa];
3. heating under vacuum with the jacket at `T_H` [K], which drives the
   remaining air out through the vent;
4. desorption with the jacket at `T_H`: steam at `T_H` and `P_L` enters the
   bottom at `steam_flow` [m³/s] and the CO2 product leaves the top through
   the vacuum pump. With `steam_flow = 0` the bed is only heated and evacuated
   (the paper's cycle B).

The `durations` [s] are those of the paper's four steps. All gas leaves the
top:

                                         ┌──V_product── air, 1 bar
    air ─────V_feed──┐                   ├──V_vent───── air, vacuum pump
                     bottom  Bed  top ───┤
    steam ──V_steam──┘                   └──V_extract── CO2 product, vacuum pump

The paper gives no permeability, only that the pressure drop stays below
1 mbar. By default the `permeability` [m²] is Blake–Kozeny's for the 7.5 mm
particles packed at the bed density, with 9% voids between them, which gives
0.05 mbar.
Axial dispersion is neglected, as in the paper, unless `dispersion` [m²/s] is
given (`nothing` for Mocca's correlation).

Heat is worth no work at the 20 °C of the air (`ambient_temperature` of the
[`Process`](@ref)), for the equivalent work of the jacket.
"""
function dac_tvsa_process(; ncells = 10, air_flow = 50e-6, steam_flow = 25e-6, P_L = 0.05e5, T_H = 368.15,
        durations = (13772.0, 30.0, 704.0, 30000.0), ramp_rate = 0.5,
        t_press = 30.0, timestepping = TimeStepping(max_dt = 100.0, first_dt = 0.01), name = "S-TVSA direct air capture", dispersion = 0.0,
        permeability = nothing,
        initial_state = Dict{Symbol, Any}(), adsorbent = load_adsorbent("apdes_nfc_stampi_bombelli_2020"),
        column = stampi_bombelli_2020_column(), gas = stampi_bombelli_2020_gas())
    feed = Feed(gas; composition = (CO2 = 0.0004, H2O = 0.0115, N2 = 0.9881), temperature = 293.15, pressure = 1e5)
    fs = Flowsheet(gas, feed)
    # The voids between the particles, at a particle density of 61 kg/m³
    ε_particles = 1 - 55.4 / 61
    add_bed!(fs, :Bed, adsorbent, column; ncells = ncells, dispersion = dispersion,
        permeability = something(permeability, compute_permeability(ε_particles, adsorbent.particle_diameter)))
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_steam, :light_inlet; at = (:Bed, :bottom),
        boundary = (pressure = P_L, temperature = T_H, composition = (H2O = 1.0,)))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_vent, :waste; at = (:Bed, :top))
    add_device!(fs, :V_extract, :extract; at = (:Bed, :top))
    t_ads, t_bd, t_heat, t_des = durations
    desorption = steam_flow > 0 ?
        Purge(t_des; flow_rate = steam_flow, jacket_temperature = T_H, name = "desorption") :
        Evacuate(t_des; to = P_L, through = :top, device = :V_extract, jacket_temperature = T_H, name = "desorption")
    cycle = Cycle([
        Pressurise(t_press),
        Adsorb(t_ads; flow_rate = air_flow),
        Blowdown(t_bd; to = P_L, device = :V_vent),
        Evacuate(t_heat; to = P_L, through = :top, device = :V_vent, jacket_temperature = T_H, name = "heating"),
        desorption,
    ])
    return Process(fs, cycle; name = name, timestepping = timestepping, ramp_rate = ramp_rate, initial_state = initial_state,
        ambient_pressure = 1e5, ambient_temperature = 293.15, heat_capacity_ratio = 1.4)
end
