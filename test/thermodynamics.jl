using Test
using Mocca
using Jutul
using StaticArrays

const ForwardDiff = Jutul.ForwardDiff

# NIST Shomate coefficients of CO2 and N2
const _SHOMATE_CO2_N2 = [(24.99735, 55.18696, -33.69137, 7.948387, -0.136638), (28.98641, 1.853978, -9.647459, 16.63537, 0.000117)]

# A model that gives only its heat capacities, so it gets the default enthalpies
struct _OnlyHeatCapacities{M} <: Mocca.AbstractGasThermodynamics
    model::M
end
Mocca.component_heat_capacities(m::_OnlyHeatCapacities, T) = Mocca.component_heat_capacities(m.model, T)

@testset "Molar mixing rule" begin
    # The mole average of the molar heat capacities, Σ yᵢ Mᵢ cₚᵢ, with the
    # constants of the components at 25 °C
    gas = Gas(["CO2", "CH4"])
    @test gas.thermodynamics isa IdealGasConstantCp
    cp = molar_heat_capacity(gas.thermodynamics, [0.4, 0.6], 298.15)
    @test cp ≈ 0.4 * 44.01e-3 * 846.0 + 0.6 * 16.043e-3 * 2226.0
    @test cp ≈ 36.3 atol = 0.05
    # Not the mass average times the mean molar mass, which gives 45.6
    @test abs(sum([0.4, 0.6] .* gas.heat_capacities) * sum([0.4, 0.6] .* gas.molar_masses) - 45.6) < 0.05
    # Constant heat capacities: the enthalpy grows linearly from the reference
    th = gas.thermodynamics
    @test molar_enthalpy(th, [0.4, 0.6], 350.0) ≈ cp * (350.0 - Mocca.THERMODYNAMIC_REFERENCE_TEMPERATURE)
    @test molar_concentration(th, 1e5, 300.0) ≈ 1e5 / (Mocca.GAS_CONSTANT * 300.0)
    @test heat_capacity_ratio(th, [0.4, 0.6], 300.0) ≈ cp / (cp - Mocca.GAS_CONSTANT)
end

@testset "Heat capacities that change with temperature" begin
    m = IdealGasShomate(_SHOMATE_CO2_N2)
    @test isbits(m) && Mocca.number_of_components(m) == 2
    @test Mocca.Physics.shomate_coefficients(m) == _SHOMATE_CO2_N2
    # NIST: CO2 37.12 and 41.34 J/(mol·K), N2 29.12 and 29.25 at 298.15 and 400 K
    @test component_heat_capacities(m, 298.15) ≈ [37.12, 29.12] atol = 0.02
    @test component_heat_capacities(m, 400.0) ≈ [41.34, 29.25] atol = 0.02
    # Water's set is fitted from 500 K; below, it stays with NIST's vapour at
    # 0.001 bar (33.612, 33.883 and 34.264 J/(mol·K) at 300, 350 and 400 K)
    w = Gas(["H2O"]; thermodynamics = :shomate).thermodynamics
    @test [component_heat_capacities(w, T)[1] for T in (300.0, 350.0, 400.0)] ≈ [33.612, 33.883, 34.264] rtol = 5e-4
    # Zero at the reference, and its derivative is the heat capacity
    @test all(iszero, component_enthalpies(m, Mocca.THERMODYNAMIC_REFERENCE_TEMPERATURE))
    y = SVector(0.15, 0.85)
    for T in (280.0, 350.0, 420.0)
        @test ForwardDiff.derivative(t -> molar_enthalpy(m, y, t), T) ≈ molar_heat_capacity(m, y, T)
    end
    # NIST: H° − H°(298.15 K) of CO2 at 400 K is 4.00 kJ/mol
    @test component_enthalpies(m, 400.0)[1] ≈ 4000.0 rtol = 3e-3
    # The default enthalpy, by quadrature of the heat capacities, agrees (the
    # E/t² term of the Shomate equation is not a polynomial)
    q = _OnlyHeatCapacities(m)
    @test Mocca.number_of_components(q) == 2
    for T in (250.0, 298.15, 360.0, 500.0)
        @test component_enthalpies(q, T) ≈ component_enthalpies(m, T) rtol = 1e-7 atol = 1e-9
    end
    # The keyword form is the same model
    @test IdealGasShomate(; A = [c[1] for c in _SHOMATE_CO2_N2], B = [c[2] for c in _SHOMATE_CO2_N2], C = [c[3] for c in _SHOMATE_CO2_N2],
        D = [c[4] for c in _SHOMATE_CO2_N2], E = [c[5] for c in _SHOMATE_CO2_N2]) == m
    @test_throws ArgumentError IdealGasShomate([(1.0, 2.0, 3.0)])

    # Automatic differentiation without allocating, as in the residuals
    cp(model, T) = molar_heat_capacity(model, y, T)
    h(model, T) = molar_enthalpy(model, y, T)
    T = ForwardDiff.Dual(330.0, 1.0)
    for model in (m, q, IdealGasConstantCp([37.0, 29.0]))
        cp(model, T), h(model, T)
        @test (@allocated cp(model, T)) == 0
        @test (@allocated h(model, T)) == 0
    end
end

@testset "Gas thermodynamics" begin
    # The constant heat capacities of the keywords as before
    gas = Gas(["CO2", "N2"]; heat_capacities = (N2 = 1000.0,))
    @test gas.thermodynamics.cp ≈ [44.01e-3 * 846.0, 28.0134e-3 * 1000.0]
    # Shomate heat capacities of the components, or given
    s = Gas(["CO2", "N2"]; thermodynamics = :shomate)
    @test s.thermodynamics == IdealGasShomate(_SHOMATE_CO2_N2)
    @test s.heat_capacities ≈ component_heat_capacities(s.thermodynamics, 298.15) ./ s.molar_masses
    s2 = Gas(["CO2", "N2"]; shomate = (N2 = (29.0, 0, 0, 0, 0),))
    @test component_heat_capacities(s2.thermodynamics, 500.0)[2] ≈ 29.0
    @test component_heat_capacities(s2.thermodynamics, 500.0)[1] ≈ component_heat_capacities(s.thermodynamics, 500.0)[1]
    # Or a model
    m = IdealGasConstantCp([30.0, 29.0])
    @test Gas(["CO2", "N2"]; thermodynamics = m).thermodynamics === m
    @test Gas(["CO2", "N2"]; thermodynamics = m).heat_capacities ≈ [30.0 / 44.01e-3, 29.0 / 28.0134e-3]
    @test occursin("IdealGasShomate", sprint(show, MIME"text/plain"(), s))
    @test_throws MoccaInputError Gas(["CO2", "N2"]; thermodynamics = m, heat_capacities = [800.0, 1000.0])
    @test_throws MoccaInputError Gas(["CO2", "N2"]; thermodynamics = IdealGasConstantCp([30.0, 29.0, 28.0]))
    @test_throws MoccaInputError Gas(["CO2", "N2"]; thermodynamics = :shomate, heat_capacities = [800.0, 1000.0])
    @test_throws MoccaInputError Gas(["CO2", "N2"]; thermodynamics = :nasa)
    @test_throws MoccaInputError Gas(["CO2", "N2"]; shomate = (N2 = (29.0, 0, 0),))
    register_component!(Component("Kr"; molar_mass = 83.798e-3, heat_capacity = 248.0, viscosity = 2.5e-5))
    @test_throws MoccaInputError Gas(["CO2", "Kr"]; thermodynamics = :shomate)
    @test Gas(["CO2", "Kr"]; shomate = (Kr = (20.786, 0, 0, 0, 0),)).thermodynamics isa IdealGasShomate

    # The adsorbed phase has the heat capacities of the gas unless given
    ads = load_adsorbent("zeolite_13x_haghpanah_2013")
    @test Mocca.Interface.adsorbed_thermodynamics(ads, s) isa IdealGasConstantCp
    free = Adsorbent(; isotherm = ads.isotherm, density = ads.density, heat_capacity = ads.heat_capacity,
        particle_diameter = ads.particle_diameter, particle_porosity = ads.particle_porosity)
    @test Mocca.Interface.adsorbed_thermodynamics(free, s) === s.thermodynamics
end

@testset "Heat capacity ratio of the pumps" begin
    # Set by the presets as in their papers, by default that of the feed
    p = four_stage_vsa_process(ncells = 4)
    @test Mocca.pump_heat_capacity_ratio(p) == 1.4
    fs = p.flowsheet
    q = Process(fs, p.cycle)
    cp = molar_heat_capacity(fs.gas.thermodynamics, fs.feed.composition, fs.feed.temperature)
    @test Mocca.pump_heat_capacity_ratio(q) ≈ cp / (cp - Mocca.GAS_CONSTANT)
    @test 1.3 < Mocca.pump_heat_capacity_ratio(q) < 1.4
    @test_throws MoccaInputError Process(fs, p.cycle; heat_capacity_ratio = 0.9)
    # Written only when set
    @test haskey(Mocca.process_to_dict(p)["options"], "heat_capacity_ratio")
    @test !haskey(Mocca.process_to_dict(q)["options"], "heat_capacity_ratio")
    @test read_process(write_process(joinpath(mktempdir(), "q.json"), q)).options.heat_capacity_ratio === nothing
end

@testset "Shomate gas in a process file" begin
    gas = Gas(["CO2", "N2"]; molar_masses = (N2 = 28e-3,), shomate = (N2 = (29.0, 0.5, 0, 0, 0),))
    p = four_stage_vsa_process(ncells = 4, gas = gas)
    d = Mocca.process_to_dict(p)
    @test d["gas"]["thermodynamics"] == "shomate" && !haskey(d["gas"], "heat_capacities")
    q = read_process(write_process(joinpath(mktempdir(), "s.json"), p))
    @test q.flowsheet.gas.thermodynamics == gas.thermodynamics
    @test q.flowsheet[:Bed].system.thermodynamics == gas.thermodynamics
    # A model of another kind cannot be written
    r = four_stage_vsa_process(ncells = 4, gas = Gas(["CO2", "N2"]; thermodynamics = _OnlyHeatCapacities(gas.thermodynamics)))
    @test_throws MoccaInputError Mocca.process_to_dict(r)
end

@testset "Heat capacities of a bed" begin
    model = Mocca.fixed_bed_model(Gas(["CO2", "N2"]; thermodynamics = :shomate), Feed(Gas(["CO2", "N2"]); composition = [0.15, 0.85]),
        load_adsorbent("zeolite_13x_haghpanah_2013"), Mocca.haghpanah_2013_column(); ncells = 3)
    sys = model.system
    @test sys.thermodynamics isa IdealGasShomate && sys.adsorbed_thermodynamics isa IdealGasConstantCp
    state0 = Mocca.setup_process_state(model; Pressure = 1e5, Temperature = [300.0, 330.0, 360.0], WallTemperature = 300.0, y = [0.5, 0.5])
    sim = Jutul.Simulator(model; state0 = state0, parameters = Jutul.setup_parameters(model))
    Jutul.update_secondary_variables!(sim.storage, model)
    s = sim.storage.state
    for cell in 1:3
        T = s.Temperature[cell]
        y = s.y[:, cell]
        q = s.AdsorbedConcentration[:, cell]
        # The gas from its composition and temperature
        @test Jutul.value(s.C_pg[cell]) ≈ molar_heat_capacity(sys.thermodynamics, Jutul.value.(y), Jutul.value(T))
        # The adsorbed phase weighted by the loadings, which are mostly CO2
        @test Jutul.value(s.C_pa[cell]) ≈ sum(Jutul.value.(q) .* sys.adsorbed_thermodynamics.cp) / sum(Jutul.value.(q))
    end
    @test Jutul.value(s.C_pg[1]) < Jutul.value(s.C_pg[3])
end

@testset "Gas at the bed's temperature leaves it there" begin
    # Argon fed at the bed's temperature into an adiabatic bed of hydrogen,
    # neither taken up by the zeolite: their molar heat capacities differ (20.8
    # and 28.8 J/(mol·K)), but nothing should change the temperature. Mocca
    # 0.2.0 moved it by 2 K, as streams brought c_p·T of their own gas and the
    # gas terms took that of the bed's.
    for thermodynamics in (nothing, :shomate)
        gas = Gas(["CO2", "Ar", "H2"]; thermodynamics = thermodynamics)
        feed = Feed(gas; composition = (Ar = 1.0,), pressure = 1e5, temperature = 298.15, light_product = (H2 = 1.0,))
        fs = Flowsheet(gas, feed)
        ads = load_adsorbent("zeolite_13x_haghpanah_2013")
        ads = Adsorbent(; isotherm = ads.isotherm, density = ads.density, heat_capacity = ads.heat_capacity,
            particle_diameter = ads.particle_diameter, particle_porosity = ads.particle_porosity)
        add_bed!(fs, :Bed, ads, Mocca.haghpanah_2013_column(); ncells = 10, thermal = Mocca.Adiabatic())
        add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
        add_device!(fs, :V_product, :product; at = (:Bed, :top))
        p = Process(fs, Cycle([Adsorb(60.0; superficial_velocity = 0.1)]); timestepping = TimeStepping(max_dt = 5.0),
            initial_state = Dict(:Bed => (pressure = 1e5, temperature = 298.15)))
        r = simulate(p)
        @test r.states[end][:Bed][:y][2, end] > 0.99
        @test maximum(maximum(abs.(s[:Bed][:Temperature] .- 298.15)) for s in r.states) < 1e-6
    end
end

@testset "Cubic heat capacities" begin
    # The coefficients of the MATLAB thermodynamics library: within 2% of the
    # NIST Shomate heat capacities from 300 to 500 K
    p = Gas(["CO2", "N2", "O2", "H2O", "CH4", "H2", "Ar"]; thermodynamics = :polynomial).thermodynamics
    s = Gas(["CO2", "N2", "O2", "H2O", "CH4", "H2", "Ar"]; thermodynamics = :shomate).thermodynamics
    @test p isa IdealGasPolynomial && isbits(p)
    for T in (300.0, 400.0, 500.0)
        @test component_heat_capacities(p, T) ≈ component_heat_capacities(s, T) rtol = 0.021
    end
    @test component_heat_capacities(p, 298.15)[1] ≈ 22.26 + 5.981e-2 * 298.15 - 3.501e-5 * 298.15^2 + 7.469e-9 * 298.15^3
    # The enthalpy is zero at the reference and its derivative is the heat capacity
    @test all(iszero, component_enthalpies(p, Mocca.THERMODYNAMIC_REFERENCE_TEMPERATURE))
    y = fill(1 / 7, 7)
    for T in (280.0, 350.0, 420.0)
        @test ForwardDiff.derivative(t -> molar_enthalpy(p, y, t), T) ≈ molar_heat_capacity(p, y, T)
    end
    @test Mocca.Physics.polynomial_coefficients(p)[2] == (28.90, -1.571e-3, 8.081e-6, -2.873e-9)
    @test IdealGasPolynomial(Mocca.Physics.polynomial_coefficients(p)) == p
    # No polynomial for CO and He
    @test_throws MoccaInputError Gas(["CO2", "CO"]; thermodynamics = :polynomial)
    @test_throws MoccaInputError Gas(["CO2", "N2"]; thermodynamics = :polynomial, heat_capacities = [800.0, 1000.0])
    # In a process file
    q = four_stage_vsa_process(ncells = 4, gas = Gas(["CO2", "N2"]; thermodynamics = :polynomial))
    @test Mocca.process_to_dict(q)["gas"]["thermodynamics"] == "polynomial"
    @test read_process(write_process(joinpath(mktempdir(), "p.json"), q)).flowsheet.gas.thermodynamics == q.flowsheet.gas.thermodynamics
end
