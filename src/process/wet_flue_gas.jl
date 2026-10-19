# VSA processes for CO2 capture from wet flue gas, after Krishnamurthy,
# Haghpanah, Rajendran and Farooq, "Simulation and Optimization of a
# Dual-Adsorbent, Two-Bed Vacuum Swing Adsorption Process for CO2 Capture from
# Wet Flue Gas", Ind. Eng. Chem. Res. 53 (2014) 14462, doi:10.1021/ie5024723.
#
# Column, particle and gas properties are from Table S2 of the paper's
# Supporting Information. It does not give the particle size; that and the
# pressure ramp rate are those of Haghpanah et al. (2013), whose model the paper
# extends.

"Components of the wet flue gas, in the order used for compositions."
const WET_FLUE_GAS_COMPONENTS = ["CO2", "N2", "H2O"]

"Feed composition of the wet flue gas: 15% CO2, 82% N2 and 3% H2O."
const WET_FLUE_GAS_FEED = [0.15, 0.82, 0.03]

"""
    wet_flue_gas_isotherm(adsorbent)

Extended dual-site Langmuir isotherm for CO2, N2 and H2O on `adsorbent`,
`:zeolite_13x` or `:silica_gel`, with the single-component parameters of Table 1
of Krishnamurthy et al. (2014). This is the paper's case 1 (competitive,
"perfect positive" form), which it uses for both adsorbents in the two-bed
process. Case 2, with CO2 parameters fitted to CO2/H2O mixture data, is a
different isotherm form and is not included.
"""
function wet_flue_gas_isotherm(adsorbent::Symbol; T0 = 298.15)
    kJ = 1e3
    if adsorbent == :zeolite_13x
        return DualSiteLangmuir(
            qsb = [3489.44, 6613.55, 10468.60], b0 = [8.65e-7, 2.50e-6, 2.35e-7], ΔUb = [-36.64, -15.82, -55.72] .* kJ,
            qsd = [2872.35, 0.0, 6434.50], d0 = [2.63e-8, 0.0, 7.99e-8], ΔUd = [-35.70, 0.0, -45.48] .* kJ,
            T0 = T0)
    elseif adsorbent == :silica_gel
        return DualSiteLangmuir(
            qsb = [2692.30, 918.90, 718.13], b0 = [3.2e-7, 1.2e-5, 5.6e-13], ΔUb = [-27.12, -12.12, -89.20] .* kJ,
            qsd = [0.0, 0.0, 52632.60], d0 = [0.0, 0.0, 7.74e-8], ΔUd = [0.0, 0.0, -39.90] .* kJ,
            T0 = T0)
    else
        error("Unknown adsorbent $adsorbent; use :zeolite_13x or :silica_gel")
    end
end

"""
Molecular diffusivities [m²/s] used for the macropore diffusion of CO2, N2 and
H2O: the CO2–N2 value of Table S2 for CO2 and N2, and the N2–H2O value for
water, which diffuses through a gas that is mostly nitrogen.
"""
const WET_FLUE_GAS_DIFFUSIVITY = SVector(1.6e-5, 1.6e-5, 2.6e-5)

"""
    wet_flue_gas_constants(adsorbent; L = 1.0, v0 = 0.45, kwarg...)

Column, particle and gas properties of Table S2 of Krishnamurthy et al. (2014)
for a bed of `adsorbent` (`:zeolite_13x` or `:silica_gel`, which differ only
in particle density) of length `L`. `v0` is the interstitial feed velocity,
which sets the axial dispersion. Other keywords override fields.
"""
function wet_flue_gas_constants(adsorbent::Symbol; L = 1.0, v0 = 0.45, kwarg...)
    ρ_s = adsorbent == :zeolite_13x ? 1130.0 : adsorbent == :silica_gel ? 1220.0 :
        error("Unknown adsorbent $adsorbent; use :zeolite_13x or :silica_gel")
    return HaghpanahConstants{Float64}(;
        L = L, r_in = 0.15, r_out = 0.16,
        Φ = 0.37, ϵ_p = 0.35, τ = 3.0, D_m = 1.6e-5,
        V0_inter = v0,
        ρ_s = ρ_s, C_ps = 1070.0,
        C_pg = SVector(1030.0, 1030.0), C_pa = SVector(1030.0, 1030.0),
        fluid_viscosity = 1.7e-5,
        K_z = 0.09, K_w = 16.0,
        ρ_w = 7800.0, C_pw = 500.0,
        h_in = 8.6, h_out = 2.5,
        T0 = 298.15, T_a = 298.15, T_feed = 298.15,
        kwarg...)
end

"""
    wet_flue_gas_system(adsorbent; constants = wet_flue_gas_constants(adsorbent),
        D_m = WET_FLUE_GAS_DIFFUSIVITY, kwarg...)

A three-component (CO2, N2, H2O) fixed bed of `adsorbent` (`:zeolite_13x` or
`:silica_gel`). Mass transfer is the macropore-diffusion linear driving force,
as the paper assumes for all three components, with a diffusivity `D_m` per
component. The heat of sorption includes every component, since water
releases most of it.
"""
function wet_flue_gas_system(adsorbent::Symbol; constants = wet_flue_gas_constants(adsorbent),
        D_m = WET_FLUE_GAS_DIFFUSIVITY, kwarg...)
    # Table S2 gives one heat capacity for the gas mixture; the adsorbed phase
    # uses the same value
    C_p = fill(constants.C_pg[1], 3)
    return FixedBed(;
        isotherm = wet_flue_gas_isotherm(adsorbent; T0 = constants.T0),
        mass_transfer = LinearDrivingForce(SVector{3}(D_m), constants.τ, constants.ϵ_p, constants.d_p),
        molecular_masses = [constants.molecularMassOfCO2, constants.molecularMassOfN2, 18.015e-3],
        component_names = WET_FLUE_GAS_COMPONENTS,
        heat_capacity_gas = C_p,
        heat_capacity_adsorbed = C_p,
        sorption_heat_all_components = true,
        kwarg...)
end

"""
    wet_flue_gas_light_product(; trace = 1e-10)

Composition of the light product used for light product pressurisation:
nitrogen with trace CO2 and H2O. The paper stores the effluent of the previous
adsorption step in a buffer; this uses its limit of pure nitrogen instead.
"""
wet_flue_gas_light_product(; trace = 1e-10) = [trace, 1 - 2trace, trace]

"""
    wet_flue_gas_initial_state(bed, pressure; temperature = 298.15)

Bed state saturated with the light product at `pressure`.
"""
function wet_flue_gas_initial_state(bed, pressure; temperature = 298.15)
    return setup_process_state(bed;
        Pressure = pressure,
        Temperature = temperature,
        WallTemperature = temperature,
        y = wet_flue_gas_light_product())
end

# Volumetric feed flow for an interstitial inlet velocity v0 [m/s]
function _interstitial_feed_rate(bed, v0)
    A = compute_column_face_area(bed, NamedTuple(Jutul.setup_parameters(bed)))
    ε = first(bed.data_domain[:porosity])
    return ε * v0 * A
end

"""
    lpp_vsa_flowsheet(; ncells = 30, v0 = 0.45,
        constants = wet_flue_gas_constants(:zeolite_13x; v0 = v0),
        system = wet_flue_gas_system(:zeolite_13x; constants = constants),
        t_lpp = 20.0, t_ads = 78.4, t_bd = 37.9, t_evac = 111.6,
        PH = 101325.0, PI = 0.08e5, PL = 0.03e5, λ = constants.λ)

The single-bed, four-step VSA cycle with light product pressurisation (LPP) on
zeolite 13X, for wet flue gas (Krishnamurthy et al. 2014, section 3). The
defaults are the operating conditions of section 3.1.

    feed ──V_feed──┐                    ┌──V_product── raffinate
                   bottom   Bed   top ──┤
    CO2 ◄─V_vacuum─┘                    └──V_lpp────── light product

| Stage          | Device    | Setting                                            |
|----------------|-----------|----------------------------------------------------|
| LPP            | V_lpp     | valve, light product ramps PL → PH, into the top   |
| adsorption     | V_feed    | flow controller at interstitial velocity `v0`      |
|                | V_product | valve, raffinate at PH                             |
| blowdown       | V_product | valve, PH → PI out of the top (co-current)         |
| evacuation     | V_vacuum  | valve, PI → PL out of the bottom (counter-current) |

The feed flow is `ε·v0·A`: the paper gives the interstitial velocity, and Mocca
works with superficial velocities.

In the paper the LPP step lasts until the bed reaches PH, at most `t_ads`, and
uses the stored effluent of the previous adsorption step. Here it lasts a
fixed `t_lpp` and uses nitrogen (see [`wet_flue_gas_light_product`](@ref)).

Returns `(fs, stages)`.
"""
function lpp_vsa_flowsheet(;
        ncells = 30,
        v0 = 0.45,
        constants = wet_flue_gas_constants(:zeolite_13x; v0 = v0),
        system = wet_flue_gas_system(:zeolite_13x; constants = constants),
        y_feed = WET_FLUE_GAS_FEED,
        T_feed = constants.T_feed,
        t_lpp = 20.0, t_ads = 78.4, t_bd = 37.9, t_evac = 111.6,
        PH = 101325.0, PI = 0.08e5, PL = 0.03e5,
        λ = constants.λ,
    )
    bed = setup_process_model(system, constants; ncells = ncells)
    names = component_names(system)
    light = wet_flue_gas_light_product()
    fs = Flowsheet()
    add_unit!(fs, :Bed, bed)
    for d in (:V_feed, :V_product, :V_lpp, :V_vacuum)
        add_unit!(fs, d, setup_device_model(names))
    end
    connect!(fs, (:Bed, :bottom) => (:V_feed, :outlet))
    connect!(fs, (:Bed, :top) => (:V_product, :inlet))
    connect!(fs, (:Bed, :top) => (:V_lpp, :outlet))
    connect!(fs, (:Bed, :bottom) => (:V_vacuum, :inlet))
    set_boundary!(fs, :V_feed, :inlet, PortCondition(pressure = PH, temperature = T_feed, composition = y_feed))
    set_boundary!(fs, :V_product, :outlet, PortCondition(pressure = PH, temperature = T_feed, composition = light))
    set_boundary!(fs, :V_lpp, :inlet, PortCondition(pressure = PH, temperature = T_feed, composition = light))
    set_boundary!(fs, :V_vacuum, :outlet, PortCondition(pressure = PL, temperature = T_feed, composition = y_feed))

    nc = Jutul.number_of_cells(bed.domain)
    C_bottom, C_top = port_conductance(bed, 1), port_conductance(bed, nc)
    ramp(p0, p1, y) = PortCondition(pressure = ExponentialRamp(p0, p1, λ), temperature = T_feed, composition = y)
    stages = [
        Stage("LPP", t_lpp; V_lpp = (law = LinearValve(C_top), inlet = ramp(PL, PH, light))),
        Stage("adsorption", t_ads;
            V_feed = (law = VolumetricFlow(_interstitial_feed_rate(bed, v0)),),
            V_product = (law = LinearValve(C_top),)),
        Stage("blowdown", t_bd; V_product = (law = LinearValve(C_top), outlet = ramp(PH, PI, light))),
        Stage("evacuation", t_evac; V_vacuum = (law = LinearValve(C_bottom), outlet = ramp(PI, PL, y_feed))),
    ]
    return (fs, stages)
end

"""
    dual_adsorbent_vsa_flowsheet(; ncells = 30,
        silica_gel_length = 0.41, zeolite_length = 1.0, v0 = 0.70,
        t_press = 20.0, t_ads = 46.20,
        silica_gel = (t_bd = 42.21, t_evac = 46.02, PI = 0.48e5, PL = 0.30e5),
        zeolite = (t_bd = 56.30, t_evac = 101.20, PI = 0.07e5, PL = 0.03e5),
        PH = 101325.0, λ = 0.5)

The dual-adsorbent, two-bed, four-step VSA process of Krishnamurthy et al.
(2014, section 4): a silica gel bed that takes up the water, in series with a
zeolite 13X bed that concentrates the CO2. The defaults are the operating
conditions for minimum energy in section 4.3.

    feed ──V_feed──┐                      ┌──V_link──┐                   ┌──V_product── raffinate
                   bottom  SilicaGel  top ┘          bottom  Zeolite  top┤
    H2O ◄─V_waste──┘                                 │                   └──V_lpp────── light product
                                          CO2 ◄─V_vacuum┘

The beds are only coupled during adsorption, when the silica gel effluent feeds
the 13X bed. Their blowdown and evacuation times and pressures differ, so the
cycle is built with [`parallel_stages`](@ref); the silica gel bed stays closed
after its evacuation until the 13X bed finishes.

| Step           | Silica gel bed                                  | 13X bed                                         |
|----------------|-------------------------------------------------|-------------------------------------------------|
| pressurisation | V_feed valve, feed ramps PL → PH                | V_lpp valve, light product ramps PL → PH        |
| adsorption     | V_feed flow controller at `v0`, V_link valve    | V_product valve, raffinate at PH                |
| blowdown       | V_waste valve, PH → PI out of the bottom        | V_product valve, PH → PI out of the top         |
| evacuation     | V_waste valve, PI → PL out of the bottom        | V_vacuum valve, PI → PL out of the bottom       |

Both beds use the properties of [`wet_flue_gas_constants`](@ref); pass
`silica_gel_constants` or `zeolite_constants` to change them. The feed flow is
`ε·v0·A`, as in [`lpp_vsa_flowsheet`](@ref).

Silica gel blows down counter-current so that the water leaves through the
feed end; 13X blows down co-current to push nitrogen out before evacuation.
Both pressurisation steps last `t_press`. The paper fixes the silica gel one
at 20 s and runs the 13X light product pressurisation until the bed reaches
PH, using the stored adsorption effluent; here it uses nitrogen (see
[`wet_flue_gas_light_product`](@ref)).

Returns `(fs, stages)`.
"""
function dual_adsorbent_vsa_flowsheet(;
        ncells = 30,
        silica_gel_length = 0.41,
        zeolite_length = 1.0,
        v0 = 0.70,
        silica_gel_constants = wet_flue_gas_constants(:silica_gel; L = silica_gel_length, v0 = v0),
        zeolite_constants = wet_flue_gas_constants(:zeolite_13x; L = zeolite_length, v0 = v0),
        silica_gel_system = wet_flue_gas_system(:silica_gel; constants = silica_gel_constants),
        zeolite_system = wet_flue_gas_system(:zeolite_13x; constants = zeolite_constants),
        y_feed = WET_FLUE_GAS_FEED,
        T_feed = 298.15,
        t_press = 20.0,
        t_ads = 46.20,
        silica_gel = (t_bd = 42.21, t_evac = 46.02, PI = 0.48e5, PL = 0.30e5),
        zeolite = (t_bd = 56.30, t_evac = 101.20, PI = 0.07e5, PL = 0.03e5),
        PH = 101325.0,
        λ = zeolite_constants.λ,
    )
    sg = setup_process_model(silica_gel_system, silica_gel_constants; ncells = ncells)
    x = setup_process_model(zeolite_system, zeolite_constants; ncells = ncells)
    names = component_names(zeolite_system)
    light = wet_flue_gas_light_product()

    fs = Flowsheet()
    add_unit!(fs, :SilicaGel, sg)
    add_unit!(fs, :Zeolite, x)
    for d in (:V_feed, :V_waste, :V_link, :V_product, :V_lpp, :V_vacuum)
        add_unit!(fs, d, setup_device_model(names))
    end
    connect!(fs, (:SilicaGel, :bottom) => (:V_feed, :outlet))
    connect!(fs, (:SilicaGel, :bottom) => (:V_waste, :inlet))
    connect!(fs, (:SilicaGel, :top) => (:V_link, :inlet))
    connect!(fs, (:Zeolite, :bottom) => (:V_link, :outlet))
    connect!(fs, (:Zeolite, :top) => (:V_product, :inlet))
    connect!(fs, (:Zeolite, :top) => (:V_lpp, :outlet))
    connect!(fs, (:Zeolite, :bottom) => (:V_vacuum, :inlet))
    set_boundary!(fs, :V_feed, :inlet, PortCondition(pressure = PH, temperature = T_feed, composition = y_feed))
    set_boundary!(fs, :V_waste, :outlet, PortCondition(pressure = silica_gel.PL, temperature = T_feed, composition = y_feed))
    set_boundary!(fs, :V_product, :outlet, PortCondition(pressure = PH, temperature = T_feed, composition = light))
    set_boundary!(fs, :V_lpp, :inlet, PortCondition(pressure = PH, temperature = T_feed, composition = light))
    set_boundary!(fs, :V_vacuum, :outlet, PortCondition(pressure = zeolite.PL, temperature = T_feed, composition = y_feed))

    nc_sg, nc_x = Jutul.number_of_cells(sg.domain), Jutul.number_of_cells(x.domain)
    C_sg_bottom, C_sg_top = port_conductance(sg, 1), port_conductance(sg, nc_sg)
    C_x_bottom, C_x_top = port_conductance(x, 1), port_conductance(x, nc_x)
    # The link joins the two half cells in series
    C_link = 1 / (1 / C_sg_top + 1 / C_x_bottom)
    ramp(p0, p1, y) = PortCondition(pressure = ExponentialRamp(p0, p1, λ), temperature = T_feed, composition = y)

    sg_stages = [
        Stage("pressurisation", t_press; V_feed = (law = LinearValve(C_sg_bottom), inlet = ramp(silica_gel.PL, PH, y_feed))),
        Stage("adsorption", t_ads;
            V_feed = (law = VolumetricFlow(_interstitial_feed_rate(sg, v0)),),
            V_link = (law = LinearValve(C_link),)),
        Stage("blowdown", silica_gel.t_bd; V_waste = (law = LinearValve(C_sg_bottom), outlet = ramp(PH, silica_gel.PI, y_feed))),
        Stage("evacuation", silica_gel.t_evac; V_waste = (law = LinearValve(C_sg_bottom), outlet = ramp(silica_gel.PI, silica_gel.PL, y_feed))),
    ]
    x_stages = [
        Stage("LPP", t_press; V_lpp = (law = LinearValve(C_x_top), inlet = ramp(zeolite.PL, PH, light))),
        Stage("adsorption", t_ads; V_product = (law = LinearValve(C_x_top),)),
        Stage("blowdown", zeolite.t_bd; V_product = (law = LinearValve(C_x_top), outlet = ramp(PH, zeolite.PI, light))),
        Stage("evacuation", zeolite.t_evac; V_vacuum = (law = LinearValve(C_x_bottom), outlet = ramp(zeolite.PI, zeolite.PL, y_feed))),
    ]
    stages = parallel_stages("silica gel" => sg_stages, "13X" => x_stages)
    return (fs, stages)
end
