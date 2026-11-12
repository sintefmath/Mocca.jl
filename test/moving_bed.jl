using Test
using Mocca
using Jutul
using StaticArrays

const R_GAS = Mocca.GAS_CONSTANT

# The activated carbon of Mondino et al. (2022), from the library
const AC = Mocca.load_adsorbent("activated_carbon_mondino_2022")
const GAS = Mocca.Gas(["CO2", "N2"]; heat_capacities = (CO2 = 846.0, N2 = 1040.0), viscosity = 1.8e-5)

@testset "Virial isotherm" begin
    iso = Mocca.isotherm(AC, GAS)
    @test iso isa Mocca.Virial
    ρ = AC.density
    # Loadings of Mondino et al. (2022, Fig. 6): feed end of the adsorption
    # section, CO2 extraction point and bottom of the cooling section
    for (T, P, y, q_expected) in (
            (303.15, 101.5e3, [0.11, 0.89], [0.35, 0.24]),
            (457.15, 90e3, [0.97, 0.03], [0.07, 0.0]),
            (291.15, 101.3e3, [0.01, 0.99], [0.06, 0.40]))
        C = SVector{2}(y .* P ./ (R_GAS * T))
        q = Mocca.compute_equilibrium(iso, C, T) ./ ρ
        @test isapprox(q, q_expected, atol = 0.01)
        @test Mocca.equilibrium_loading(AC, GAS; pressure = P, temperature = T, composition = y) ≈ q
        # The loading satisfies the implicit isotherm equation
        A = iso.A0 .+ iso.A1 ./ T
        B = iso.B0 .+ iso.B1 ./ T
        S = sum(q)
        g = (A .* S .+ sum(A .* q)) ./ 2 .+ (B .* S^2 .+ 2 * S * sum(B .* q)) ./ 3
        K = iso.K∞ .* exp.(iso.ΔH ./ (R_GAS * T))
        @test isapprox(q ./ K .* exp.(g), C .* R_GAS .* T, rtol = 1e-10)
    end
    # Henry's law at low pressure
    T = 320.0
    C = SVector(1e-6, 1e-6)
    K = iso.K∞ .* exp.(iso.ΔH ./ (R_GAS * T))
    @test isapprox(Mocca.compute_equilibrium(iso, C, T), ρ .* K .* C .* R_GAS .* T, rtol = 1e-6)
    @test Mocca.compute_enthalpy(iso, C, T) ≈ -[30.006e3, 21.934e3]

    # Derivatives through the Newton solve match finite differences
    f(x) = Mocca.compute_equilibrium(iso, SVector(x[1], x[2]), x[3])
    x = [4.4, 35.8, 303.15]
    J = Jutul.ForwardDiff.jacobian(f, x)
    for k in 1:3
        h = zeros(3); h[k] = 1e-6 * x[k]
        @test isapprox(J[:, k], (f(x .+ h) .- f(x .- h)) ./ (2h[k]), rtol = 1e-6)
    end

    # The Langmuir isotherm of the MATLAB parameter set holds far less CO2 at the feed
    matlab = Mocca.Adsorbent(; isotherm = Mocca.langmuir("CO2" => (qs = 0.4, b0 = 2.7e-7, ΔU = -27e3),
            "N2" => (qs = 0.02, b0 = 5.5e-7, ΔU = -18e3); basis = :mass),
        density = AC.density, heat_capacity = AC.heat_capacity, particle_diameter = AC.particle_diameter, particle_porosity = AC.particle_porosity)
    q = Mocca.equilibrium_loading(matlab, GAS; pressure = 101.5e3, temperature = 303.15, composition = [0.11, 0.89])
    @test isapprox(q[1], 0.020, atol = 1e-3)
end

@testset "Virial isotherm data" begin
    data = AC.isotherm
    @test data isa Mocca.VirialData
    # Aliases, and omitted terms are zero
    v = Mocca.virial("CO2" => (Kinf = 1e-9, dH = 30e3, A0 = 1.0), "N2" => (K∞ = 2e-9, ΔH = 20e3))
    @test v.K∞ == [1e-9, 2e-9] && v.ΔH == [30e3, 20e3] && v.A0 == [1.0, 0.0] && v.B1 == [0.0, 0.0]
    @test_throws Mocca.MoccaInputError Mocca.virial("CO2" => (K∞ = 1e-9, ΔH = 30e3, C0 = 1.0))
    # Every component of the gas needs parameters, and loadings per kg need the density
    @test_throws Mocca.MoccaInputError Mocca.build_isotherm(Mocca.virial("CO2" => (K∞ = 1e-9, ΔH = 30e3)), GAS; density = 900.0)
    @test_throws Mocca.MoccaInputError Mocca.build_isotherm(data, GAS)
    # In the order of the gas
    swapped = Mocca.isotherm(AC, Mocca.Gas(["N2", "CO2"]))
    @test swapped.K∞ == reverse(Mocca.isotherm(AC, GAS).K∞)
    # Written and read back
    path = Mocca.save_adsorbent(joinpath(mktempdir(), "ac.json"), AC)
    back = Mocca.load_adsorbent(path)
    for k in fieldnames(Mocca.VirialData)
        @test getfield(back.isotherm, k) == getfield(data, k)
    end
end

@testset "Macropores of the activated carbon" begin
    # Knudsen diffusion in the 0.525 µm macropores, in the order of the gas
    @test AC.macropore_radius == 5.25e-7
    mt = Mocca.mass_transfer(AC, GAS)
    @test Mocca.has_knudsen(mt) && mt.molar_masses == GAS.molar_masses
    swapped = Mocca.mass_transfer(AC, Mocca.Gas(["N2", "CO2"]))
    @test Mocca.pore_diffusivity(swapped, 303.15) ≈ reverse(Mocca.pore_diffusivity(mt, 303.15))
    @test all(Mocca.pore_diffusivity(mt, 303.15) .< 1.6e-5)
    # Written and read back, also as an explicit mass transfer
    back = Mocca.load_adsorbent(Mocca.save_adsorbent(joinpath(mktempdir(), "ac.json"), AC))
    @test back.macropore_radius == AC.macropore_radius
    explicit = Mocca.Adsorbent("explicit"; isotherm = AC.isotherm, density = AC.density, heat_capacity = AC.heat_capacity,
        particle_diameter = AC.particle_diameter, particle_porosity = AC.particle_porosity, mass_transfer = mt)
    back = Mocca.load_adsorbent(Mocca.save_adsorbent(joinpath(mktempdir(), "explicit.json"), explicit))
    @test Mocca.pore_diffusivity(back.mass_transfer, 400.0) == Mocca.pore_diffusivity(mt, 400.0)
    # The radius sets the default mass transfer, so not both
    @test_throws Mocca.MoccaInputError Mocca.Adsorbent(; isotherm = AC.isotherm, density = 904.0, heat_capacity = 880.0,
        particle_diameter = 7e-4, particle_porosity = 0.5, macropore_radius = 5e-7, mass_transfer = mt)
    @test_throws Mocca.MoccaInputError Mocca.Adsorbent(; isotherm = AC.isotherm, density = 904.0, heat_capacity = 880.0,
        particle_diameter = 7e-4, particle_porosity = 0.5, macropore_radius = 0.0)
end

@testset "Fixed molar flow and back pressure controller" begin
    st = (InletPressure = [1.1e5], OutletPressure = [1.5e5], InletTemperature = [300.0], OutletTemperature = [350.0])
    @test Mocca.device_flow(Mocca.FixedMolarFlow(3.0), st, 0.0) == 3.0
    # Draws from the inlet whatever the outlet pressure
    bpc = Mocca.BackPressureController(1e5, 1e-3)
    @test Mocca.device_flow(bpc, st, 0.0) ≈ 1e-3 * 1e4 * 1.1e5 / (R_GAS * 300.0)
    # Reverses below the setpoint, with the outlet as upstream
    st_low = merge(st, (InletPressure = [0.9e5],))
    @test Mocca.device_flow(bpc, st_low, 0.0) ≈ -1e-3 * 1e4 * 1.5e5 / (R_GAS * 350.0)
    model = Mocca.setup_device_model(["CO2", "N2"])
    for law in (Mocca.FixedMolarFlow(3.0), bpc)
        v = zeros(Jutul.vectorization_length(law, model, :law, :all))
        Jutul.vectorize_force!(v, model, law, :law, :all)
        @test Jutul.devectorize_force(law, model, v, nothing, :law, :all) == law
    end
end

# A small lab-scale moving bed of the activated carbon (Virial isotherm)
# The gas in the particles' pores moves with the sorbent
function lab_moving_bed(; thermal = Mocca.Adiabatic(), ncells = 6, length = 0.3, radius = 0.025,
        sorbent_mass_flow = 0.01, heat_exchanger = nothing, particle_porosity = AC.particle_porosity,
        limiter = Mocca.VanAlbada(), solid_inlet_temperature = 303.15, kwarg...)
    sys = Mocca.MovingBed(; isotherm = Mocca.isotherm(AC, GAS), mass_transfer = Mocca.mass_transfer(AC, GAS),
        molecular_masses = GAS.molar_masses, component_names = GAS.names, heat_capacity_gas = GAS.heat_capacities,
        heat_capacity_adsorbed = GAS.heat_capacities, thermal = thermal, limiter = limiter)
    return Mocca.setup_moving_bed_model(sys;
        ncells = ncells, length = length, radius = radius, porosity = 0.5,
        particle_diameter = AC.particle_diameter, adsorbent_density = AC.density,
        adsorbent_heat_capacity = AC.heat_capacity, sorbent_mass_flow = sorbent_mass_flow,
        solid_inlet_temperature = solid_inlet_temperature, dispersion = 1e-4, thermal_conductivity = 0.2,
        particle_porosity = particle_porosity, heat_exchanger = heat_exchanger, kwarg...)
end

lab_tubes(; T_in, ṁ = 0.02) = Mocca.tube_bundle(area_density = 229.0, sorbent_side_htc = 100.0,
    fluid_mass_flow = ṁ, fluid_inlet_temperature = T_in)

const Y_FEED = [0.11, 0.89]
const FEED_RATE = 1e-4 / 0.02976  # mol/s, about 5 L/min

# Feed into the bottom of `unit`, gas out of the top to atmospheric pressure
function add_feed_and_vent!(fs, unit)
    names = GAS.names
    Mocca.add_unit!(fs, :V_feed, Mocca.setup_device_model(names))
    Mocca.add_unit!(fs, :V_vent, Mocca.setup_device_model(names))
    Mocca.connect!(fs, (unit, :bottom) => (:V_feed, :outlet))
    Mocca.connect!(fs, (unit, :top) => (:V_vent, :inlet))
    Mocca.set_boundary!(fs, :V_feed, :inlet, Mocca.PortCondition(pressure = 1.02e5, temperature = 303.15, composition = Y_FEED))
    Mocca.set_boundary!(fs, :V_vent, :outlet, Mocca.PortCondition(pressure = 1.01325e5, temperature = 303.15, composition = [1e-6, 1 - 1e-6]))
    return (V_feed = (law = Mocca.FixedMolarFlow(FEED_RATE),),
        V_vent = (law = Mocca.LinearValve(Mocca.port_conductance(fs[unit], Jutul.number_of_cells(fs[unit].domain))),))
end

function simulate_flowsheet(fs, settings, init; duration = 20.0, max_dt = 2.0, parameters = Mocca.setup_flowsheet_parameters(fs))
    mm = Mocca.setup_flowsheet_model(fs)
    state0 = Mocca.setup_flowsheet_state(fs; init...)
    forces, dt = Mocca.setup_schedule(fs, mm, [Mocca.Stage("run", duration, settings...)]; max_dt = max_dt, first_dt = 0.1)
    # States, lengths and forces of every step the solver took
    states, ts, step_forces = Mocca.Operation._simulate_steps(mm, state0, parameters, forces, dt; nonlinear_tolerance = 1e-10)
    return (; mm, state0, parameters, forces, dt, states, ts, step_forces)
end

@testset "Limiters" begin
    limiters = (Mocca.VanAlbada(), Mocca.VanLeer(), Mocca.Minmod())
    @test Mocca.limited_slope(Mocca.Upwind(), 1.0, 2.0) == 0
    for l in limiters
        # Second order on a straight line, and first order at extrema
        @test Mocca.limited_slope(l, 0.3, 0.3) ≈ 0.3
        @test Mocca.limited_slope(l, -0.3, -0.3) ≈ -0.3
        @test Mocca.limited_slope(l, 1.0, -2.0) == 0 && Mocca.limited_slope(l, -1.0, 0.5) == 0
        @test Mocca.limited_slope(l, 0.0, 0.0) == 0 && Mocca.limited_slope(l, 0.0, 1.0) == 0
        # Within the TVD region: the sign of a, and at most 2·min(|a|, |b|)
        for (a, b) in ((1.0, 0.1), (1.0, 10.0), (0.2, 0.7), (-2.0, -0.5))
            s = Mocca.limited_slope(l, a, b)
            @test sign(s) == sign(a) && abs(s) <= 2 * min(abs(a), abs(b)) + 1e-14
        end
        # Smooth in the slopes, as Newton and the adjoints need, where they agree in sign
        g = Jutul.ForwardDiff.gradient(x -> Mocca.limited_slope(l, x[1], x[2]), [0.4, 0.9])
        @test all(isfinite, g)
        @test all(isfinite, Jutul.ForwardDiff.gradient(x -> Mocca.limited_slope(l, x[1], x[2]), [0.0, 0.0]))
    end
    # A linear profile is reconstructed exactly at the faces, away from the ends
    φ(c) = 2.0 + 0.5 * c
    @test Mocca.Balances._sorbent_face(Mocca.VanAlbada(), φ, 3, 6) ≈ (φ(3) + φ(2)) / 2
    @test Mocca.Balances._sorbent_face(Mocca.Upwind(), φ, 3, 6) == φ(3)
    # First order next to the inlet and at the outlet
    @test Mocca.Balances._sorbent_face(Mocca.VanAlbada(), φ, 6, 6) == φ(6)
    @test Mocca.Balances._sorbent_face(Mocca.VanAlbada(), φ, 1, 6) == φ(1)
end

@testset "A limiter sharpens a thermal front" begin
    # Hot sorbent moving down into a bed at 30 °C: the front is sharper with
    # van Albada than with first-order upwinding, and the temperatures stay
    # between those of the bed and the inlet
    # The sorbent and its pores enter in equilibrium with nitrogen at 80 °C
    T_in = 353.15
    C_in = [1e-6, 1 - 1e-6] .* 1.01325e5 ./ (R_GAS * T_in)
    q_in = collect(Mocca.compute_equilibrium(Mocca.isotherm(AC, GAS), C_in, T_in))
    function front(limiter)
        bed = lab_moving_bed(ncells = 60, limiter = limiter, solid_inlet_temperature = T_in, sorbent_mass_flow = 0.02,
            solid_inlet_loading = q_in, solid_inlet_gas_concentration = C_in)
        fs = Mocca.UnitNetwork()
        Mocca.add_unit!(fs, :Bed, bed)
        settings = add_feed_and_vent!(fs, :Bed)
        settings = merge(settings, (V_feed = (law = Mocca.Closed(),),))
        init = (Bed = Mocca.setup_unit_state(bed; pressure = 1.01325e5, temperature = 303.15, composition = [1e-6, 1 - 1e-6]),)
        r = simulate_flowsheet(fs, pairs(settings), init; duration = 6.0, max_dt = 0.25)
        T = r.states[end][:Bed][:Temperature]
        θ = (T .- 303.15) ./ 50
        width = count(x -> 0.1 < x < 0.9, θ)
        return (; width, T)
    end
    up, va = front(Mocca.Upwind()), front(Mocca.VanAlbada())
    @test up.width > 0 && va.width < 0.75 * up.width
    @test minimum(va.T) > 303.15 - 0.5 && maximum(va.T) < 353.15 + 0.5
end

@testset "A moving bed without sorbent flow is a fixed bed" begin
    wall = (r_out = 0.03, density = 7800.0, heat_capacity = 500.0, conductivity = 16.0,
        inner_htc = 10.0, outer_htc = 2.0, ambient_temperature = 303.15)
    # A fixed bed has no gas in particle pores
    moving = lab_moving_bed(sorbent_mass_flow = 0.0, wall = wall, particle_porosity = 0.0)
    sys = moving.system
    fixed_sys = Mocca.FixedBed(; isotherm = sys.isotherm, mass_transfer = sys.mass_transfer,
        molecular_masses = sys.molecular_masses, component_names = sys.component_names,
        thermodynamics = sys.thermodynamics, adsorbed_thermodynamics = sys.adsorbed_thermodynamics,
        thermal = Mocca.Adiabatic())
    fixed = Mocca.setup_process_model(fixed_sys, moving.data_domain)
    init = (pressure = 1.01325e5, temperature = 303.15, composition = [1e-6, 1 - 1e-6])
    results = map((moving, fixed)) do bed
        fs = Mocca.UnitNetwork()
        Mocca.add_unit!(fs, :Bed, bed)
        settings = add_feed_and_vent!(fs, :Bed)
        simulate_flowsheet(fs, pairs(settings), (Bed = Mocca.setup_unit_state(bed; init...),); duration = 30.0)
    end
    a, b = results[1].states[end][:Bed], results[2].states[end][:Bed]
    # CO2 has entered the bed and is being taken up
    @test maximum(a[:y][1, :]) > 0.05
    for k in (:Pressure, :Temperature, :AdsorbedConcentration, :y)
        @test isapprox(a[k], b[k], rtol = 1e-6)
    end
end

# Adsorption, preheating, regeneration and cooling sections with the sorbent
# going round. The water that cools the sorbent in the cooler goes on to
# preheat it, and the regenerator is heated by hot water.
const LOOP_SECTIONS = (:Adsorber, :Preheater, :Regenerator, :Cooler)

function four_section_loop(; inlet_loading = zeros(2), inlet_temperature = 303.15, inlet_gas = zeros(2), preheater_fluid_inlet = 300.0)
    fs = Mocca.UnitNetwork()
    inlet = (solid_inlet_loading = inlet_loading, solid_inlet_temperature = inlet_temperature, solid_inlet_gas_concentration = inlet_gas)
    hx_section(T_in, ṁ) = lab_moving_bed(; thermal = Mocca.HeatExchanger(), ncells = 5, length = 0.112, radius = 0.063,
        heat_exchanger = lab_tubes(T_in = T_in, ṁ = ṁ), inlet...)
    Mocca.add_unit!(fs, :Adsorber, lab_moving_bed(; inlet...))
    # The recovered water carries about as much heat per kelvin as the sorbent
    Mocca.add_unit!(fs, :Preheater, hx_section(preheater_fluid_inlet, 0.0025))
    Mocca.add_unit!(fs, :Regenerator, hx_section(430.0, 0.02))
    Mocca.add_unit!(fs, :Cooler, hx_section(288.15, 0.0025))
    for (a, b) in zip(LOOP_SECTIONS, (LOOP_SECTIONS[2:end]..., LOOP_SECTIONS[1]))
        Mocca.connect_solids!(fs, a => b)
    end
    Mocca.connect_heat_transfer_fluid!(fs, :Cooler => :Preheater)
    settings = Dict(pairs(add_feed_and_vent!(fs, :Adsorber)))
    # Gas released in the preheater is vented, the regenerator gives the
    # product, and the cooler takes in nitrogen as the sorbent cools
    names = GAS.names
    for d in (:V_preheat, :V_product, :V_cool)
        Mocca.add_unit!(fs, d, Mocca.setup_device_model(names))
    end
    Mocca.connect!(fs, (:Preheater, :top) => (:V_preheat, :inlet))
    Mocca.connect!(fs, (:Regenerator, :bottom) => (:V_product, :inlet))
    Mocca.connect!(fs, (:Cooler, :bottom) => (:V_cool, :outlet))
    Mocca.set_boundary!(fs, :V_preheat, :outlet, Mocca.PortCondition(pressure = 1.01325e5, temperature = 303.15, composition = Y_FEED))
    Mocca.set_boundary!(fs, :V_product, :outlet, Mocca.PortCondition(pressure = 0.95e5, temperature = 303.15, composition = [0.9, 0.1]))
    Mocca.set_boundary!(fs, :V_cool, :inlet, Mocca.PortCondition(pressure = 1.01325e5, temperature = 303.15, composition = [1e-6, 1 - 1e-6]))
    nc = Jutul.number_of_cells(fs[:Preheater].domain)
    settings[:V_preheat] = (law = Mocca.LinearValve(Mocca.port_conductance(fs[:Preheater], nc)),)
    settings[:V_product] = (law = Mocca.LinearValve(Mocca.port_conductance(fs[:Regenerator], 1)),)
    settings[:V_cool] = (law = Mocca.LinearValve(Mocca.port_conductance(fs[:Cooler], 1)),)
    init = (
        Adsorber = Mocca.setup_unit_state(fs[:Adsorber]; pressure = 1.02e5, temperature = 303.15, composition = Y_FEED),
        Preheater = Mocca.setup_unit_state(fs[:Preheater]; pressure = 1.01325e5, temperature = 330.0, composition = [0.3, 0.7]),
        Regenerator = Mocca.setup_unit_state(fs[:Regenerator]; pressure = 0.95e5, temperature = 400.0, composition = [0.8, 0.2]),
        Cooler = Mocca.setup_unit_state(fs[:Cooler]; pressure = 1.01325e5, temperature = 340.0, composition = [0.1, 0.9]),
    )
    return (fs, settings, init)
end

@testset "Moving bed flowsheet validation" begin
    fs, = four_section_loop()
    @test_throws ErrorException Mocca.connect_solids!(fs, :Adsorber => :Cooler)       # the cooler is already fed
    @test_throws ErrorException Mocca.connect_solids!(fs, :Adsorber => :Adsorber)
    @test_throws ErrorException Mocca.connect_solids!(fs, :V_feed => :Adsorber)       # not a moving bed
    @test_throws ErrorException Mocca.connect_heat_transfer_fluid!(fs, :Adsorber => :Cooler)  # no heat exchanger
    @test_throws ErrorException lab_moving_bed(thermal = Mocca.HeatExchanger())     # tubes missing
    # Units passing sorbent to each other must agree on its flow
    net = Mocca.UnitNetwork()
    Mocca.add_unit!(net, :A, lab_moving_bed(sorbent_mass_flow = 0.01))
    Mocca.add_unit!(net, :B, lab_moving_bed(sorbent_mass_flow = 0.02))
    @test_throws ErrorException Mocca.connect_solids!(net, :A => :B)
    # Sorbent and fluid connections name their own ports, apart from the gas ports
    c = Mocca.SorbentConnection(:Cooler, :Adsorber)
    @test Mocca.connection_ports(c) == (Mocca.PortRef(:Cooler, :sorbent_outlet), Mocca.PortRef(:Adsorber, :sorbent_inlet))
    fs, = four_section_loop()
    mm = Mocca.setup_flowsheet_model(fs)
    # Gas: two cross terms per connection plus one for the energy balance;
    # sorbent: loading, pore gas and energy per connection; fluid: one
    @test length(mm.cross_terms) == 5 * 3 + 4 * 3 + 1
end

@testset "Moving bed loop conserves each component" begin
    # Long enough for the sorbent to go round several times (about 220 s)
    fs, settings, init = four_section_loop()
    r = simulate_flowsheet(fs, pairs(settings), init; duration = 900.0, max_dt = 30.0)
    inventory(state) = sum(Mocca.component_inventory(fs[u], state[u], r.parameters[u]) for u in LOOP_SECTIONS)
    totals(d) = Mocca.stream_totals(r.states, r.ts, d)
    n_in = totals(:V_feed) .+ totals(:V_cool)
    n_out = totals(:V_vent) .+ totals(:V_preheat) .+ totals(:V_product)
    change = inventory(r.states[end]) .- inventory(r.state0)
    @test isapprox(change, n_in .- n_out, rtol = 1e-6, atol = 1e-9)
    s = r.states[end]
    # The regenerator gives a product enriched in CO2
    n_product = totals(:V_product)
    @test n_product[1] / sum(n_product) > 3 * Y_FEED[1]
    # Water warmed in the cooler preheats the sorbent leaving the adsorber
    T_recovered = s[:Cooler][:HXFluidTemperature][end]
    @test T_recovered > 288.15 + 5
    @test s[:Preheater][:HXFluidTemperature][1] < T_recovered
    @test s[:Preheater][:Temperature][1] > s[:Adsorber][:Temperature][1]
end

@testset "Connected inlets ignore their inlet parameters" begin
    fs_a, settings, init = four_section_loop()
    fs_b, = four_section_loop(inlet_loading = [500.0, 100.0], inlet_temperature = 450.0, inlet_gas = [30.0, 5.0], preheater_fluid_inlet = 400.0)
    a = simulate_flowsheet(fs_a, pairs(settings), init; duration = 10.0)
    b = simulate_flowsheet(fs_b, pairs(settings), init; duration = 10.0)
    for u in LOOP_SECTIONS, k in keys(fs_a[u].primary_variables)
        @test a.states[end][u][k] ≈ b.states[end][u][k]
    end
end

@testset "Moving bed adjoint gradients" begin
    # Gradients of the CO2 product with respect to parameters of the moving
    # beds and of the streams between them match central finite differences
    fs, settings, init = four_section_loop()
    prm = Mocca.setup_flowsheet_parameters(fs)
    r = simulate_flowsheet(fs, pairs(settings), init; duration = 10.0, max_dt = 1.0, parameters = prm)
    # The CO2 product, plus temperatures that respond within the short run to
    # the water entering the cooler (through the fluid stream to the
    # preheater) and to the adsorber's heat capacity
    product = Mocca.StreamObjective(:V_product; scale = 1e6)
    G(model, state, dt, step_info, forces) = product(model, state, dt, step_info, forces) +
        dt * (state[:Preheater][:HXFluidTemperature][1] + state[:Adsorber][:Temperature][1])
    function objective(p)
        states, ts, sf = Mocca.Operation._simulate_steps(r.mm, r.state0, p, r.forces, r.dt; nonlinear_tolerance = 1e-10)
        return Mocca.Operation._evaluate_objective(G, r.mm, states, ts, sf)
    end
    grad = Jutul.solve_adjoint_sensitivities(r.mm, r.states, r.ts, G;
        state0 = r.state0, parameters = prm, forces = r.step_forces, info_level = -1)
    sections = LOOP_SECTIONS
    # Each with a relative step large enough to stand above the solver's
    # tolerance in the finite differences
    perturbations = [
        # Sorbent flow, the same in every section
        ([(u, :SorbentMassFlow) for u in sections], 1e-4),
        # Heat exchange, and the water entering the cooler, which then preheats
        ([(:Regenerator, :SorbentSideHeatTransferCoeff)], 1e-4),
        ([(:Regenerator, :HXFluidInletTemperature)], 1e-4),
        ([(:Cooler, :HXFluidInletTemperature)], 2.5e-3),
        ([(:Preheater, :HXFluidMassFlow), (:Cooler, :HXFluidMassFlow)], 1e-4),
        ([(:Regenerator, :TubeOuterArea)], 1e-4),
        ([(:Adsorber, :AdsorbentHeatCapacity)], 1e-4),
        # The gas carried in the pores of the sorbent
        ([(u, :ParticlePorosity) for u in sections], 1e-4),
    ]
    for (group, rel) in perturbations
        pp, pm = deepcopy(prm), deepcopy(prm)
        dG_adj = 0.0
        for (u, k) in group
            h = rel .* abs.(prm[u][k])
            pp[u][k] .+= h
            pm[u][k] .-= h
            dG_adj += sum(grad[u][k] .* h)
        end
        dG_fd = (objective(pp) - objective(pm)) / 2
        @test abs(dG_fd) > 0
        @test isapprox(dG_adj, dG_fd, rtol = 1e-5)
    end
    # Inlet parameters of connected units have no effect
    for u in sections
        @test all(iszero, grad[u][:SolidInletTemperature])
        @test all(iszero, grad[u][:SolidInletLoading])
        @test all(iszero, grad[u][:SolidInletGasConcentration])
    end
    @test all(iszero, grad[:Preheater][:HXFluidInletTemperature])
end

# The MBTSA process of Mondino et al. (2022) on a coarse grid, as in
# examples/mbtsa_mondino_2022.jl. Without `paper_dispersion` the sections get
# Mocca's dispersion instead of the paper's correlation.
function mbtsa_process(; cells_per_metre = 40, duration = 2400.0, gas = GAS, paper_dispersion = true)
    y_feed = [0.11, 0.89]
    feed = Mocca.Feed(gas; composition = y_feed, temperature = 303.15, pressure = 101.5e3,
        light_product = Mocca.light_product_composition(y_feed; trace = 1e-6))
    M(y) = sum(y .* gas.molar_masses)
    sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
    heights = (0.9, 0.4, 0.8, 0.4, 0.6)
    diameters = (12.8, 13.3, 13.3, 13.3, 13.3)
    porosity = (0.7, 0.5, 0.5, 0.5, 0.5)
    packing = (0.05, 0.0, 0.0, 0.0, 0.0)
    ṁ_s = 650.0
    # Dispersion and conductivity of the paper's correlations
    area(i) = π * (diameters[i] / 2)^2
    F_feed = 55.9 / M(y_feed)
    velocity = (F_feed * R_GAS * 303.15 / 101.5e3 / area(1), 0.05, 0.05, 0.05, 0.05) ./ porosity
    Re(i) = 1.2 * velocity[i] * porosity[i] * AC.particle_diameter / 1.8e-5
    Sc = 1.8e-5 / (1.2 * 1.6e-5)
    function hx(i, ṁ, T_in)
        v = ṁ_s / (AC.density * (1 - porosity[i] - packing[i]) * area(i))
        return Mocca.tube_bundle(area_density = 229.0, sorbent_side_htc = Mocca.moving_bed_sorbent_htc(v),
            fluid_side_htc = 10000.0, fluid_mass_flow = ṁ, fluid_inlet_temperature = T_in)
    end
    hxs = (nothing, hx(2, 160.6, 303.15), hx(3, 267.6, 460.15), hx(4, 160.6, 303.15), hx(5, 267.6, 283.15))
    fs = Mocca.Flowsheet(gas, feed)
    for (i, name) in enumerate(sections)
        column = Mocca.Column(length = heights[i], radius = diameters[i] / 2, porosity = porosity[i],
            axial_conductivity = 0.026 * (7 + 0.5 * 0.7 * Re(i)))
        Mocca.add_moving_bed!(fs, name, AC, column; sorbent_mass_flow = ṁ_s, ncells = max(4, round(Int, cells_per_metre * heights[i])),
            packing_fraction = packing[i], heat_exchanger = hxs[i],
            dispersion = paper_dispersion ? 1.6e-5 / porosity[i] * (20 + 0.5 * Sc * Re(i)) : nothing)
    end
    for (i, name) in enumerate(sections)
        Mocca.connect_solids!(fs, name => sections[mod1(i + 1, 5)])
    end
    Mocca.connect_heat_transfer_fluid!(fs, :Precooling => :Preheating)
    Mocca.connect_heat_transfer_fluid!(fs, :Preheating => :Precooling)
    Mocca.add_device!(fs, :V_feed, :feed; at = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_raffinate, :product; at = (:Adsorption, :top), boundary = (pressure = 101.325e3,))
    Mocca.add_device!(fs, :V_purge, :link; from = (:Adsorption, :top), to = (:Cooling, :bottom))
    Mocca.add_device!(fs, :V_link, :link; from = (:Cooling, :top), to = (:Precooling, :bottom))
    Mocca.add_device!(fs, :V_recycle_preheat, :link; from = (:Preheating, :top), to = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_recycle_precool, :link; from = (:Precooling, :top), to = (:Adsorption, :bottom))
    Mocca.add_device!(fs, :V_product, :extract; at = (:Desorption, :bottom), boundary = (pressure = 90e3, temperature = 460.15, composition = [0.97, 0.03]))
    C(d) = Mocca.default_conductance(fs, d)
    settings = Dict(
        :V_feed => Mocca.FixedMolarFlow(F_feed),
        :V_raffinate => Mocca.LinearValve(C(:V_raffinate)),
        :V_purge => Mocca.FixedMolarFlow(0.15 * 55.9 / M(feed.light_product)),
        :V_link => Mocca.LinearValve(C(:V_link)),
        :V_recycle_preheat => Mocca.BackPressureController(101.325e3, C(:V_recycle_preheat)),
        :V_recycle_precool => Mocca.BackPressureController(101.325e3, C(:V_recycle_precool)),
        :V_product => Mocca.LinearValve(C(:V_product)))
    initial = Dict(
        :Adsorption => (temperature = 303.15, pressure = 101.5e3, composition = y_feed),
        :Preheating => (temperature = 383.15, pressure = 101.325e3, composition = [0.5, 0.5]),
        :Desorption => (temperature = 453.15, pressure = 90e3, composition = [0.95, 0.05]),
        :Precooling => (temperature = 393.15, pressure = 101.325e3, composition = [0.5, 0.5]),
        :Cooling => (temperature = 313.15, pressure = 101.325e3, composition = [0.02, 0.98]))
    return Mocca.SteadyProcess(fs, settings; name = "MBTSA", duration = duration,
        timestepping = Mocca.TimeStepping(max_dt = 60.0, first_dt = 0.05), heat_capacity_ratio = 1.4, initial_state = initial)
end

@testset "Steady processes" begin
    p = mbtsa_process()
    fs = p.flowsheet
    @test Mocca.sorbent_path(fs) == [:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling]
    @test Mocca.sorbent_path(fs, :Desorption)[end] == :Preheating
    @test Mocca.moving_sorbent_units(fs) == Mocca.sorbent_path(fs)
    diagram = sprint(Mocca.print_flowsheet_diagram, fs)
    @test occursin("moving bed of activated carbon", diagram) && occursin("sorbent to Preheating", diagram)
    # A cyclic process runs a cycle of steps, which a moving bed does not have
    @test_throws Mocca.MoccaInputError Mocca.Process(fs, Mocca.Cycle([Mocca.Adsorb(10.0)]))
    # A steady process is a process: one stage, and the places it has and
    # none it does not
    @test p isa Mocca.Process && Mocca.SteadyResult <: Mocca.ProcessResult
    @test only(Mocca.stages(p)).settings == Mocca.operating_stage(p).settings
    q = deepcopy(p)
    Mocca.add_variable!(q, :water, (:parameter, :Desorption, :HXFluidMassFlow); bounds = (200.0, 300.0))
    @test q[:water] == 267.6
    q[:water] = 250.0
    @test all(==(250.0), q.parameters[:Desorption][:HXFluidMassFlow]) && q.parameter_settings[:Desorption][:HXFluidMassFlow] == 250.0
    @test Mocca.variable_values(q) == [250.0]
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(q, :t, ("adsorption", :duration))
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(q, :P, (:feed, :pressure))
    # Variables of a steady process are written with it
    Mocca.add_variable!(q, :extraction, (:setting, :V_product, :outlet, :pressure))
    @test q[:extraction] == 90e3
    q[:extraction] = 85e3
    @test q.settings[:V_product].outlet.pressure == 85e3 && q.settings[:V_product].outlet.temperature == 460.15
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(q, :bad, (:setting, :V_link, :inlet, :pressure))
    r = Mocca.read_process(Mocca.write_process(joinpath(mktempdir(), "steady.json"), q))
    @test Mocca.variable_names(r) == [:water, :extraction] && r[:water] == 250.0 && r[:extraction] == 85e3
    # A device setting: a field of the device's flow law
    Mocca.add_variable!(q, :purge, (:setting, :V_purge, :rate))
    Mocca.add_variable!(q, :setpoint, Mocca.SettingPlace(:V_recycle_preheat, :setpoint))
    rate = q[:purge]
    q[:purge] = 1.5 * rate
    @test q.settings[:V_purge].law == Mocca.FixedMolarFlow(1.5 * rate) && p.settings[:V_purge].law.rate == rate
    q[:setpoint] = 101.0e3
    @test q.settings[:V_recycle_preheat].law.setpoint == 101.0e3 && q.settings[:V_recycle_preheat].law.conductance == p.settings[:V_recycle_preheat].law.conductance
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(q, :bad, (:setting, :V_purge, :conductance))
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(q, :bad, (:setting, :Adsorption, :rate))
    @test_throws Mocca.MoccaInputError Mocca.add_variable!(Mocca.four_stage_vsa_process(ncells = 4), :bad, (:setting, :V_feed, :rate))
    # Settings name devices, and boundaries only unconnected ports
    @test_throws Mocca.MoccaInputError Mocca.SteadyProcess(fs, Dict(:Adsorption => Mocca.Closed()))
    @test_throws Mocca.MoccaInputError Mocca.SteadyProcess(fs, Dict(:V_link => (law = Mocca.Closed(),
        inlet = Mocca.PortCondition(pressure = 1e5, temperature = 300.0, composition = [0.5, 0.5]))))
    @test_throws Mocca.MoccaInputError Mocca.SteadyProcess(fs, Dict(:V_feed => Mocca.FixedMolarFlow(1.0)); duration = 0.0)

    # Written to a file and read back, it is the same process
    path = Mocca.write_process(joinpath(mktempdir(), "mbtsa.json"), p)
    q = Mocca.read_process(path)
    @test q isa Mocca.SteadyProcess
    @test Mocca.process_to_dict(q) == Mocca.process_to_dict(p)
    # The limiter is written by name
    d = Mocca.process_to_dict(p)
    @test all(u -> u["limiter"] == "van_albada", d["units"])
    d["units"][1]["limiter"] = "upwind"
    @test Mocca.Interface.process_from_dict(d).flowsheet[:Adsorption].system.limiter isa Mocca.Upwind
    d["units"][1]["limiter"] = "superbee"
    @test_throws Mocca.MoccaInputError Mocca.Interface.process_from_dict(d)
    a = Mocca.simulate(p; duration = 5.0)
    b = Mocca.simulate(q; duration = 5.0)
    @test sum(a.timesteps) ≈ 5.0
    for u in Mocca.sorbent_path(fs), k in (:Pressure, :Temperature, :AdsorbedConcentration, :HXFluidTemperature)
        haskey(a.states[end][u], k) || continue
        @test a.states[end][u][k] ≈ b.states[end][u][k]
    end
end

@testset "MBTSA process of Mondino et al. (2022)" begin
    p = mbtsa_process()
    r = Mocca.steady_state(p; max_spans = 6)
    @test r.converged
    perf = Mocca.performance(r)
    # Steady: what enters leaves, and the profiles have stopped changing
    @test isapprox(perf.extract .+ perf.product, perf.feed, rtol = 1e-6)
    capture(s) = Mocca.recovery(Mocca.stream_rates(s, :V_product), Mocca.stream_rates(s, :V_feed))
    @test capture(r.states[end - 5]) ≈ perf.recovery rtol = 1e-4
    # The recovery loop moves heat from the hot sorbent to the cold
    @test perf.heat_duty[:Preheating] ≈ -perf.heat_duty[:Precooling] rtol = 1e-8
    # Only the desorption water heats from outside, and only the cooling water cools
    @test perf.heating ≈ perf.heat_duty[:Desorption]
    @test perf.cooling ≈ -perf.heat_duty[:Cooling]
    # The heat put in leaves through the cooling water, less the heat of
    # sorption of the CO2 taken out and the sensible heat of the products
    @test isapprox(perf.heating, perf.cooling, rtol = 0.02)
    # Cycle time of the sorbent, 4.2 min in the paper
    @test isapprox(perf.circulation_time, 4.2 * 60, rtol = 0.02)
    @test perf.sorbent_inventory ≈ perf.circulation_time * 650.0
    # On this coarse grid; at 200 cells per metre the process gives 99.7%
    # purity, 93.1% capture and 5.53 MJ/kg, within 0.4 points of finer grids
    # (see the example)
    @test perf.purity ≈ 0.9884 atol = 1e-3
    @test perf.recovery ≈ 0.9121 atol = 1e-3
    @test perf.specific_heat_duty ≈ 5.669e6 rtol = 1e-3

    # The recycle and purge blowers raise the gas they move; the valve that
    # passes the purge on lets it down
    fs, o, state = p.flowsheet, p.options, r.states[end]
    one_step(d) = Mocca.blower_energy([state], [1.0], d)
    @test one_step(:V_recycle_preheat) > 0 && one_step(:V_recycle_precool) > 0 && one_step(:V_purge) > 0
    @test one_step(:V_link) == 0
    @test perf.blower_power ≈ sum(one_step, (:V_purge, :V_link, :V_recycle_preheat, :V_recycle_precool))
    pumps = Mocca.vacuum_pump_energy([state], [1.0], :V_product; discharge_pressure = o.ambient_pressure) +
        Mocca.compressor_energy([state], [1.0], :V_feed; suction_pressure = o.ambient_pressure)
    @test perf.power ≈ pumps + perf.blower_power
    # The equivalent work of the hot water, ṁ·c_p·[(T_in − T_out) − T₀·ln(T_in/T_out)],
    # about a third of the heat at 187 °C
    prm = p.parameters[:Desorption]
    T_in, T_out = prm[:HXFluidInletTemperature][1], state[:Desorption][:HXFluidTemperature][end]
    mc = prm[:HXFluidMassFlow][1] * prm[:HXFluidHeatCapacity][1]
    @test perf.equivalent_work ≈ mc * ((T_in - T_out) - 298.15 * log(T_in / T_out)) rtol = 1e-12
    @test perf.equivalent_work ≈ Mocca.equivalent_work(fs, state, p.parameters, :Desorption)
    @test 0.3 < perf.equivalent_work / perf.heating < 0.4
    M = GAS.molar_masses[1]
    @test perf.specific_equivalent_work ≈ (perf.power + perf.equivalent_work) / 3.6e6 / (perf.extract[1] * M / 1000)
    # A measure of the user's own: the equivalent work per tonne of CO2 with
    # the heat weighted by a second-law efficiency of 0.75
    Mocca.add_measure!(p, :weighted_work, s -> (s.energy + 0.75 * s.equivalent_work) / (s.key_extract * M))
    @test Mocca.performance(r).weighted_work ≈ (perf.power + 0.75 * perf.equivalent_work) / (perf.extract[1] * M)
end

@testset "Steady state gradients" begin
    # At steady state, by the adjoint of one long step whose fixed point the
    # steady state is, against central differences of steady states found by
    # running the process
    p = mbtsa_process(cells_per_metre = 10)
    sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
    Mocca.add_variable!(p, :sorbent_flow, ((:parameter, u, :SorbentMassFlow) for u in sections)...)
    Mocca.add_variable!(p, :T_cold, (:parameter, :Cooling, :HXFluidInletTemperature))
    Mocca.add_variable!(p, :purge, (:setting, :V_purge, :rate))
    # The height of the desorption section, and the pressure the product is
    # drawn to, a boundary condition of the flowsheet
    Mocca.add_variable!(p, :height, (:geometry, :Desorption, :length))
    Mocca.add_variable!(p, :extraction, (:setting, :V_product, :outlet, :pressure))
    ss = Mocca.steady_state(p; tol = 1e-9, max_spans = 20)
    @test ss.converged
    # The next solve starts from this steady state
    again = Mocca.steady_state(p; tol = 1e-9, max_spans = 20)
    @test again.converged && length(again.history) == 1
    measures = [:recovery, :specific_heat_duty, :specific_equivalent_work]
    g = Mocca.performance_gradient(p; measures = measures, steady_state = ss, progress = false)
    @test g.metrics.recovery == Mocca.performance(ss).recovery
    x_ss = Mocca.restart_state(p.model, ss.states[end])
    function value(n, x)
        x0 = p[n]
        p[n] = x
        r = Mocca.steady_state(p; tol = 1e-11, max_spans = 40, state0 = x_ss)
        # Measured before the variable is reset, since the heat duty reads
        # the process's parameters
        perf = Mocca.performance(r)
        p[n] = x0
        return perf
    end
    for n in Mocca.variable_names(p)
        x0 = p[n]
        h = 1e-4 * x0
        up, down = value(n, x0 + h), value(n, x0 - h)
        for m in measures
            @test g.gradients[m][n] ≈ (up[m] - down[m]) / 2h rtol = 1e-5
        end
    end
end

@testset "Moving beds with transport properties that change" begin
    # A gas whose viscosity and diffusivities follow Chapman–Enskog: each
    # section scales them by temperature and composition, and its dispersion
    # follows them by Wakao and Funazkri
    gas = Mocca.Gas(["CO2", "N2"]; heat_capacities = (CO2 = 846.0, N2 = 1040.0), transport = :chapman_enskog)
    p = mbtsa_process(cells_per_metre = 10, gas = gas, paper_dispersion = false)
    sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
    for u in sections
        sys = p.flowsheet[u].system
        @test sys.transport.model isa Mocca.ChapmanEnskogTransport
        @test sys.dispersion isa WakaoFunazkri
        @test haskey(p.parameters[u], :DispersionScale) && !haskey(p.parameters[u], :AxialDispersion)
    end
    # Given a dispersion, a section keeps it as its parameter
    q = mbtsa_process(cells_per_metre = 10, gas = gas)
    @test q.flowsheet[:Desorption].system.dispersion isa ConstantDispersion
    @test haskey(q.parameters[:Desorption], :AxialDispersion)

    ss = Mocca.steady_state(p; tol = 1e-9, max_spans = 20)
    @test ss.converged
    # The gas is more viscous in the hot desorption section than at the feed
    function viscosity_factors(u)
        st, tr = ss.states[end][u], p.flowsheet[u].system.transport
        return [Mocca.Physics.viscosity_factor(tr, st[:y][:, c], st[:Temperature][c]) for c in eachindex(st[:Temperature])]
    end
    @test all(>(1.1), viscosity_factors(:Desorption))
    @test all(x -> isapprox(x, 1, atol = 0.05), viscosity_factors(:Adsorption))
    # And the transport properties matter
    constant = Mocca.steady_state(mbtsa_process(cells_per_metre = 10, paper_dispersion = false); tol = 1e-9, max_spans = 20)
    @test !isapprox(Mocca.performance(ss).recovery, Mocca.performance(constant).recovery, rtol = 1e-4)

    # Gradients through the scaled viscosity, mass transfer and dispersion,
    # against central differences (see check_gradients)
    Mocca.add_variable!(p, :sorbent_flow, ((:parameter, u, :SorbentMassFlow) for u in sections)...)
    Mocca.add_variable!(p, :T_hot, (:parameter, :Desorption, :HXFluidInletTemperature))
    Mocca.add_variable!(p, :purge, (:setting, :V_purge, :rate))
    Mocca.add_variable!(p, :height, (:geometry, :Desorption, :length))
    Mocca.add_variable!(p, :disp, (:parameter, :Adsorption, :DispersionScale))
    checks = Mocca.check_gradients(p; measures = (:recovery, :specific_heat_duty, :specific_equivalent_work), io = nothing)
    @test all(r -> r.status == :pass, checks.rows)
end

@testset "Adjoints of the blower and heat duty integrals" begin
    Op = Mocca.Operation
    p = mbtsa_process(cells_per_metre = 10)
    fs, model, prm = p.flowsheet, p.model, p.parameters
    forces, dt = Mocca.Interface._span_schedule(p, 30.0)
    x0 = Mocca.initial_state(p)
    G = (heat = Mocca.HeatDutyObjective(fs, :Desorption), loop = Mocca.HeatDutyObjective(fs, :Preheating),
        blower = Mocca.BlowerObjective(:V_recycle_preheat), work = Mocca.EquivalentWorkObjective(fs, :Desorption),
        valve = Mocca.BlowerObjective(:V_link))
    function values(forces, prm)
        local s, t, f = Op._simulate_steps(model, x0, prm, forces, dt)
        return map(g -> Op._evaluate_objective(g, model, s, t, f, prm), G)
    end
    v = values(forces, prm)
    @test v.blower > 0 && v.valve == 0
    # Summed over the steps, the heat duty integral is the heat duty times the time
    states, ts, sf = Op._simulate_steps(model, x0, prm, forces, dt)
    @test v.heat ≈ sum(Mocca.heat_duty(fs, s, prm, :Desorption) * h for (s, h) in zip(states, ts)) rtol = 1e-12
    # Against central differences in the purge rate, a setting of the stage
    with_purge(rate) = [merge(f, Dict(:V_purge => merge(f[:V_purge], (law = Mocca.FixedMolarFlow(rate),)))) for f in forces]
    r0 = forces[1][:V_purge].law.rate
    h = 1e-4 * r0
    up, down = values(with_purge(r0 + h), prm), values(with_purge(r0 - h), prm)
    for k in keys(G)
        g = Mocca.force_gradients(model, x0, prm, forces, dt, G[k])
        @test g.objective ≈ v[k]
        @test g.forces[1][:V_purge].law.rate ≈ (up[k] - down[k]) / 2h rtol = 1e-5 atol = 1e-8
    end
    # And in the hot water flow, a parameter read by the integrals themselves
    m0 = prm[:Desorption][:HXFluidMassFlow][1]
    hm = 1e-4 * m0
    with_flow(m) = (q = deepcopy(prm); q[:Desorption][:HXFluidMassFlow] .= m; q)
    up, down = values(forces, with_flow(m0 + hm)), values(forces, with_flow(m0 - hm))
    for k in (:heat, :work, :blower)
        g = Jutul.solve_adjoint_sensitivities(model, states, ts, G[k]; state0 = x0, parameters = prm, forces = sf, info_level = -1)
        @test sum(g[:Desorption][:HXFluidMassFlow]) ≈ (up[k] - down[k]) / 2hm rtol = 1e-5
    end
end
