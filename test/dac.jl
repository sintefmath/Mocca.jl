using Test
using Mocca
using Mocca.Units
using Jutul
using StaticArrays

const M = Mocca

@testset "Toth isotherm with water on a GAB isotherm" begin
    # Saturation pressure of water: 2.339 kPa at 20 °C, 1 atm at 100 °C
    @test isapprox(water_vapour_pressure(293.15), 2339.0, rtol = 0.01)
    @test isapprox(water_vapour_pressure(373.15), 101325.0, rtol = 0.005)

    gas = Gas(["CO2", "H2O", "N2"])
    ads = load_adsorbent("apdes_nfc_stampi_bombelli_2020")
    iso = M.isotherm(ads, gas)
    @test iso isa TothGAB
    # The equations of Stampi-Bombelli et al. (2020) in mol/kg and kPa, with
    # the parameters of their Table 2
    function paper(T, p_c, p_w)
        x = p_w / (water_vapour_pressure(T) / 1e3)
        q_w = 36.48 * 0.1489 * 0.5751x / ((1 - 0.5751x) * (1 + (0.1489 - 1) * 0.5751x))
        T0 = 296.0
        ns = 2.38 / (1 - 0.0061q_w)
        b = 70.74e4 * exp(57047 / (M.GAS_CONSTANT * T0) * (T0 / T - 1)) * (1 + 28.907q_w)
        t = 0.4148 - 1.606 * (1 - T0 / T)
        return (ns * b * p_c / (1 + (b * p_c)^t)^(1 / t), q_w)
    end
    for (T, P, y) in ((293.15, 1e5, [0.0004, 0.0115, 0.9881]), (368.15, 5e3, [0.05, 0.95, 0.0]), (323.15, 1e5, [0.0004, 0.0, 0.9996]))
        q = equilibrium_loading(ads, gas; pressure = P, temperature = T, composition = y)
        q_c, q_w = paper(T, y[1] * P / 1e3, y[2] * P / 1e3)
        @test q[1] ≈ q_c rtol = 1e-10
        @test q[2] ≈ q_w rtol = 1e-10 atol = 1e-14
        @test q[3] == 0
    end
    # Humidity raises the uptake of CO2, and heat lowers it
    humid = equilibrium_loading(ads, gas; pressure = 1e5, temperature = 293.15, composition = (CO2 = 0.0004, H2O = 0.0115, N2 = 0.9881))
    dry = equilibrium_loading(ads, gas; pressure = 1e5, temperature = 293.15, composition = (CO2 = 0.0004, N2 = 0.9996))
    hot = equilibrium_loading(ads, gas; pressure = 1e5, temperature = 368.15, composition = (CO2 = 0.0004, H2O = 0.0115, N2 = 0.9881))
    @test humid[1] > dry[1] > hot[1]
    # The heats of adsorption are those given, in gas order
    @test compute_enthalpy(iso, SVector(1.0, 1.0, 1.0), 300.0) == SVector(-57000.0, -49000.0, 0.0)
    # By component name: the same loadings for a gas in another order, and none
    # for a gas without water
    other = M.isotherm(ads, Gas(["N2", "CO2", "H2O"]))
    C = SVector(10.0, 0.016, 0.47)
    @test compute_equilibrium(other, C, 300.0) ≈ compute_equilibrium(iso, SVector(C[2], C[3], C[1]), 300.0)[[3, 1, 2]]
    dry_iso = M.isotherm(ads, Gas(["CO2", "N2"]))
    @test dry_iso.water == 0
    @test compute_equilibrium(dry_iso, SVector(0.016, 40.0), 296.0)[1] ≈ compute_equilibrium(iso, SVector(0.016, 0.0, 40.0), 296.0)[1]
    # Saturation capacities per m³ of adsorbent, for the cycle change
    @test M.saturation_capacity(iso) ≈ SVector(2.38, 36.48, 0.0) .* ads.density
    # The gas may hold far less CO2 than in other beds, and only in this bed
    @test M.minimum_mole_fraction(iso) == 1e-20
    @test M.minimum_mole_fraction(M.isotherm(load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_gas())) == 1e-12
    # A negative partial pressure, which Newton's method may visit, gives a
    # finite loading
    @test all(isfinite, compute_equilibrium(iso, SVector(-1e-6, 0.4, 40.0), 293.15))
    @test_throws MoccaInputError toth_gab(ns0 = 2.38, b0 = 707.4, t0 = 0.41, ΔH0 = -57e3, reference_temperature = 296.0,
        cm = 36.48, cG = 0.15, K = 0.58, heats = (CO2 = -57e3,), water = "CO2")
    @test_throws MoccaInputError toth_gab(ns0 = -1.0, b0 = 707.4, t0 = 0.41, ΔH0 = -57e3, reference_temperature = 296.0,
        cm = 36.48, cG = 0.15, K = 0.58, heats = (CO2 = -57e3,))
end

@testset "Lumped linear driving force" begin
    gas = Gas(["CO2", "H2O", "N2"])
    ads = load_adsorbent("apdes_nfc_stampi_bombelli_2020")
    mt = M.mass_transfer(ads, gas)
    @test mt isa LumpedLinearDrivingForce
    @test mt.k == SVector(2e-4, 2e-3, 1.0)
    @test compute_mass_transfer_rate(mt, nothing, SVector(1.0, 2.0, 0.0), SVector(3.0, 1.0, 0.0)) ≈ SVector(4e-4, -2e-3, 0.0)
    @test M.build_mass_transfer(lumped_ldf("H2O" => 3.0; default = 0.5), gas).k == SVector(0.5, 3.0, 0.5)
    @test_throws MoccaInputError lumped_ldf("CO2" => 0.0)
    @test_throws MoccaInputError lumped_ldf("CO2" => 1.0, "CO2" => 2.0)
    # The adsorbent, isotherm and mass transfer, through a file
    path = joinpath(mktempdir(), "nfc.json")
    save_adsorbent(path, ads)
    again = load_adsorbent(path)
    @test M.Interface.isotherm_to_dict(again.isotherm) == M.Interface.isotherm_to_dict(ads.isotherm)
    @test again.mass_transfer.k == ads.mass_transfer.k && again.mass_transfer.default == ads.mass_transfer.default
    @test equilibrium_loading(again, gas; pressure = 1e5, temperature = 300.0, composition = [0.0004, 0.01, 0.9896]) ==
        equilibrium_loading(ads, gas; pressure = 1e5, temperature = 300.0, composition = [0.0004, 0.01, 0.9896])
end

# A small bed of 13X, heated by its jacket while it is closed between
# adsorption and evacuation
function heated_process(; T_j = 400.0, ncells = 4)
    gas = M.haghpanah_2013_gas()
    fs = Flowsheet(gas, Feed(gas; composition = [0.15, 0.85], pressure = 1bar))
    column = Column(; length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37, outer_heat_transfer = 500.0)
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), column; ncells = ncells)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    cycle = Cycle([Pressurise(5.0), Adsorb(5.0; velocity = 1.0), Idle(10.0; jacket_temperature = T_j, name = "heating"),
        Evacuate(10.0; to = 0.1bar)])
    return Process(fs, cycle; timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.5))
end

@testset "Heating jackets" begin
    p = heated_process()
    sts = M.stages(p)
    @test sts[3].settings[:Bed] == (jacket_temperature = 400.0,)
    @test all(!haskey(st.settings, :Bed) for st in sts[[1, 2, 4]])
    @test occursin("Bed jacket at 126.9 °C", sprint(show, MIME"text/plain"(), p))
    # A jacket on a step with no single bed, or a cold one, is reported
    @test_throws MoccaInputError Process(p.flowsheet, Cycle([Adsorb(5.0; velocity = 1.0), Idle(5.0; jacket_temperature = -1.0)]))

    r = simulate(p)
    t = cumsum(r.timesteps)
    Tw(k) = maximum(r.states[k][:Bed][:WallTemperature])
    k_ads, k_heat, k_evac = findlast(<=(10.0 + 1e-9), t), findlast(<=(20.0 + 1e-9), t), length(t)
    # The wall warms while the jacket is hot, and cools when it is back at
    # ambient temperature
    @test Tw(k_heat) > Tw(k_ads) + 5
    @test Tw(k_evac) < Tw(k_heat)
    # Heat flows from the jacket to the wall only while it is hot
    Q = jacket_heat(r)
    @test Q > 0
    # The jacket heat is the heating of the process, with its specific measure
    perf = performance(r)
    @test perf.heating ≈ Q rtol = 1e-12
    @test haskey(perf, :specific_heat_duty) && haskey(perf, :specific_equivalent_work)
    # Its equivalent work: the heat weighted by the Carnot factor at the jacket
    W = M.integral_values((W = JacketHeatObjective(:Bed; ambient_temperature = 298.15),), p.model, r.states, r.timesteps,
        p.parameters; forces = M.Interface._step_forces(r)).W
    @test W ≈ Q * (1 - 298.15 / 400.0) rtol = 1e-12
    # A process without a jacket has none of the heat measures
    @test !haskey(M.performance_measures(four_stage_vsa_process(ncells = 4)), :specific_heat_duty)
    prm = p.parameters[:Bed]
    Q_check = sum(r.timesteps[k] * sum(prm[:WallAreaOut] .* 500.0 .* (400.0 .- r.states[k][:Bed][:WallTemperature]))
        for k in eachindex(t) if 10.0 < t[k] - r.timesteps[k] / 2 < 20.0)
    G_ends = 2 * 16.0 * prm[:WallCrossSectionArea][1] / (prm[:CellDx][1] / 2)
    @test Q > Q_check
    @test Q < Q_check + G_ends * 10.0 * (400.0 - 298.15)
    # Without the jacket the heat is zero and the run is that of ambient
    # surroundings
    cold = simulate(heated_process(T_j = 298.15))
    @test abs(jacket_heat(cold)) < 1e-3 * Q

    # The jacket temperature of a step through a file
    path = write_process(joinpath(mktempdir(), "heated.json"), p)
    q = read_process(path)
    @test M.stages(q)[3].settings[:Bed] == (jacket_temperature = 400.0,)
end

@testset "Gradient with respect to a jacket temperature" begin
    p = heated_process()
    forces, dt = M.cycle_schedule(p)
    x0, prm = M.initial_state(p), p.parameters
    G = M.StreamObjective(:V_vacuum; scale = 1e3)
    g = M.force_gradients(p.model, x0, prm, forces, dt, G)
    @test isnothing(g.forces[2][:Bed].jacket_temperature)
    function objective(T)
        sets = unique(objectid, forces)
        new = IdDict{Any, Any}(f => (i == 3 ? merge(f, Dict(:Bed => Jutul.setup_forces(p.model[:Bed]; jacket_temperature = T))) : f)
            for (i, f) in enumerate(sets))
        f = [new[s] for s in forces]
        states, ts, sf = M.Operation._simulate_steps(p.model, x0, prm, f, dt)
        return M.Operation._evaluate_objective(G, p.model, states, ts, sf)
    end
    @test g.objective ≈ objective(400.0)
    h = 0.01
    fd = (objective(400.0 + h) - objective(400.0 - h)) / (2h)
    @test fd > 0  # a hotter bed gives up more CO2 to the vacuum
    @test isapprox(g.forces[3][:Bed].jacket_temperature, fd, rtol = 1e-5)
end

@testset "Direct air capture by S-TVSA" begin
    p = dac_tvsa_process()
    @test [st.name for st in M.stages(p)] == ["pressurisation", "adsorption", "blowdown", "heating", "desorption"]
    @test M.stages(p)[5].settings[:V_steam].law == M.VolumetricFlow(25e-6)
    # Newton's method, the default, settles it
    r = cyclic_steady_state(p)
    @test r.converged
    perf = performance(r)
    n = perf.extract
    M_CO2 = p.flowsheet.gas.molar_masses[1]
    # The production rate is within the paper's range, and the CO2 product dry
    # and pure
    @test 1.9e-5 < n[1] * M_CO2 / perf.cycle_time * 3600 < 2.6e-5
    i_N2 = findfirst(==("N2"), p.flowsheet.gas.names)
    @test n[1] / (n[1] + n[i_N2]) > 0.996
    @test perf.heating > 0
    # Without steam the bed is only heated and evacuated
    b = dac_tvsa_process(; steam_flow = 0.0)
    @test M.stages(b)[5].settings[:V_extract].law isa M.LinearValve
    @test !haskey(M.stages(b)[5].settings, :V_steam)
end
