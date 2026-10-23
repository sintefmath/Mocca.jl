using Test
using Mocca
using Mocca.Units
using Jutul
using StaticArrays

const ForwardDiff = Jutul.ForwardDiff

# A model whose properties do not change, so that a reference with other
# values gives fixed factors on every cell
struct _FixedTransport <: Mocca.AbstractGasTransport end
Mocca.mixture_viscosity(::_FixedTransport, y, T) = one(T)
Mocca.mixture_thermal_conductivity(::_FixedTransport, y, T) = one(T)
Mocca.diffusivity_scale(::_FixedTransport, T, P) = one(T / P)

@testset "Component transport properties" begin
    gas = Gas(["CO2", "N2", "H2O"]; transport = :sutherland)
    t = gas.transport
    @test t isa SutherlandTransport && isbits(t)
    # NIST Chemistry WebBook at 400 K and 0.01 bar
    @test Mocca.Physics.component_viscosities(t, 400.0) ≈ [19.621e-6, 22.200e-6, 13.354e-6] rtol = 0.025
    @test Mocca.Physics.component_thermal_conductivities(t, 400.0) ≈ [0.024674, 0.032782, 0.026435] rtol = 0.02
    # At 25 °C the mixture viscosity is that of the components mixed as before
    y = [0.15, 0.82, 0.03]
    @test mixture_viscosity(t, y, 298.15) == Mocca.gas_viscosity(Gas(["CO2", "N2", "H2O"]), y)
    # One component: its own conductivity
    @test mixture_thermal_conductivity(t, [1.0, 0.0, 0.0], 350.0) ≈ Mocca.Physics.component_thermal_conductivities(t, 350.0)[1]
    # Viscosity and conductivity rise with temperature; diffusivity as T^1.75/P
    @test mixture_viscosity(t, y, 350.0) > mixture_viscosity(t, y, 300.0)
    @test mixture_thermal_conductivity(t, y, 350.0) > mixture_thermal_conductivity(t, y, 300.0)
    @test Mocca.diffusivity_scale(t, 400.0, 0.5e5) / Mocca.diffusivity_scale(t, 300.0, 1e5) ≈ (4 / 3)^1.75 * 2
    # With Knudsen diffusion in series: molecular at high pressure, and at low
    # pressure independent of it and going as √T
    PK = 0.2e5
    k = sutherland_transport(["CO2", "N2", "H2O"]; knudsen_pressure = PK)
    D(T, P) = Mocca.diffusivity_scale(k, T, P)
    @test D(350.0, 0.5e5) ≈ 1 / (0.5e5 * (298.15 / 350)^1.75 + PK * sqrt(298.15 / 350))
    @test D(298.15, PK) ≈ 1 / 2PK   # where the two are equal
    @test D(400.0, 1.0) / D(300.0, 1.0) ≈ sqrt(4 / 3) rtol = 1e-4
    @test D(300.0, 1e9) / D(300.0, 2e9) ≈ 2 rtol = 1e-4
    @test k.μ₀ == t.μ₀ && k.S == t.S
    @test_throws MoccaInputError sutherland_transport(["CO2"]; knudsen_pressure = -1.0)
    @test_throws ArgumentError SutherlandTransport(; viscosity = [1.5e-5, 1.7e-5], sutherland = [270.0], thermal_conductivity = [0.017],
        conductivity_exponent = [1.3], molar_masses = [44e-3])

    # Automatic differentiation without allocating
    ys = SVector(0.15, 0.82, 0.03)
    T = ForwardDiff.Dual(330.0, 1.0)
    f(t, T) = (mixture_viscosity(t, ys, T), mixture_thermal_conductivity(t, ys, T), Mocca.diffusivity_scale(t, T, 1e5))
    f(t, T)
    @test (@allocated f(t, T)) == 0
    h = 1e-3
    @test ForwardDiff.partials(f(t, T)[1])[1] ≈ (mixture_viscosity(t, ys, 330.0 + h) - mixture_viscosity(t, ys, 330.0 - h)) / 2h rtol = 1e-6
end

@testset "Factors against the reference" begin
    t = Gas(["CO2", "N2"]; transport = :sutherland).transport
    y, T, P = [0.15, 0.85], 298.15, 1e5
    r = TransportReference(t, y, T, P)
    @test Mocca.viscosity_factor(r, y, T) == 1 && Mocca.conductivity_factor(r, y, T) == 1 && Mocca.diffusivity_factor(r, T, P) == 1
    @test Mocca.viscosity_factor(r, y, 350.0) ≈ mixture_viscosity(t, y, 350.0) / mixture_viscosity(t, y, T)
    @test Mocca.diffusivity_factor(r, 350.0, 0.5e5) ≈ (350 / 298.15)^1.75 * 2
    c = TransportReference(ConstantTransport())
    @test Mocca.is_constant_transport(c) && !Mocca.is_constant_transport(r)
    @test TransportReference(ConstantTransport(), y, T, P) == c
    @test Mocca.viscosity_factor(c, y, 400.0) === 1.0 && Mocca.diffusivity_factor(c, 400.0, 1e4) === 1.0

    # The linear driving force scales its molecular diffusivity
    mt = LinearDrivingForce(D_m = 1.6e-5, τ = 3.0, ϵ_p = 0.35, d_p = 2e-3)
    C, q, qs = SVector(5.0, 30.0), SVector(100.0, 50.0), SVector(800.0, 60.0)
    @test compute_mass_transfer_rate(mt, C, q, qs, 2.0) ≈
        compute_mass_transfer_rate(LinearDrivingForce(D_m = 3.2e-5, τ = 3.0, ϵ_p = 0.35, d_p = 2e-3), C, q, qs)
end

@testset "Gas transport" begin
    @test Gas(["CO2", "N2"]).transport === ConstantTransport()
    @test Gas(["CO2", "N2"]; transport = _FixedTransport()).transport === _FixedTransport()
    @test_throws MoccaInputError Gas(["CO2", "N2"]; transport = :kinetic)
    register_component!(Component("Ne"; molar_mass = 20.18e-3, heat_capacity = 1030.0, viscosity = 3.2e-5))
    @test_throws MoccaInputError Gas(["CO2", "Ne"]; transport = :sutherland)
    register_component!(Component("Ne"; molar_mass = 20.18e-3, heat_capacity = 1030.0, viscosity = 3.2e-5, sutherland = 56.0,
        thermal_conductivity = 0.049, conductivity_exponent = 0.69))
    @test Gas(["CO2", "Ne"]; transport = :sutherland).transport.S[2] == 56.0
    # The viscosity at the feed follows the feed temperature with a model
    g = Gas(["CO2", "N2"]; transport = :sutherland)
    @test Mocca.gas_viscosity(g, [0.15, 0.85], 350.0) ≈ mixture_viscosity(g.transport, [0.15, 0.85], 350.0)
    @test Mocca.gas_viscosity(Gas(["CO2", "N2"]), [0.15, 0.85], 350.0) == Mocca.gas_viscosity(Gas(["CO2", "N2"]), [0.15, 0.85])
    @test Mocca.gas_viscosity(Gas(["CO2", "N2"]; transport = :sutherland, viscosity = 1.7e-5), [0.15, 0.85], 350.0) == 1.7e-5

    # In a process file
    p = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = :sutherland))
    d = Mocca.process_to_dict(p)
    @test d["gas"]["transport"] == "sutherland"
    q = read_process(write_process(joinpath(mktempdir(), "t.json"), p))
    @test q.flowsheet.gas.transport == p.flowsheet.gas.transport
    @test q.flowsheet[:Bed].system.transport == p.flowsheet[:Bed].system.transport
    @test !haskey(Mocca.process_to_dict(four_stage_vsa_process(ncells = 4))["gas"], "transport")
    kn = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = sutherland_transport(["CO2", "N2"];
        molar_masses = [44.01e-3, 28e-3], knudsen_pressure = 0.1e5)))
    @test Mocca.process_to_dict(kn)["gas"]["transport"]["knudsen_pressure"] == 0.1e5
    @test read_process(write_process(joinpath(mktempdir(), "k.json"), kn)).flowsheet.gas.transport == kn.flowsheet.gas.transport
    r = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = _FixedTransport()))
    @test_throws MoccaInputError Mocca.process_to_dict(r)
end

@testset "A bed with transport properties that change" begin
    # The bed refers its parameters to the feed, and has the factors only
    # when its transport changes
    p = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = :sutherland))
    bed = p.flowsheet[:Bed]
    feed = p.flowsheet.feed
    @test bed.system.transport == TransportReference(p.flowsheet.gas.transport, feed.composition, feed.temperature, feed.pressure)
    @test haskey(bed.secondary_variables, :ViscosityFactor) && haskey(bed.secondary_variables, :ConductivityFactor)
    @test bed.secondary_variables[:AdsorptionMassTransfer] isa Mocca.ScaledAdsorptionMassTransfer
    plain = four_stage_vsa_process(ncells = 4).flowsheet[:Bed]
    @test !haskey(plain.secondary_variables, :ViscosityFactor) && !haskey(plain.secondary_variables, :ConductivityFactor)
    @test plain.secondary_variables[:AdsorptionMassTransfer] isa Mocca.AdsorptionMassTransfer
    # Without an energy balance there is no conductivity to scale
    iso = Mocca.fixed_bed_model(p.flowsheet.gas, feed, load_adsorbent("zeolite_13x_haghpanah_2013"), Mocca.haghpanah_2013_column();
        ncells = 3, thermal = Mocca.Isothermal(300.0))
    @test haskey(iso.secondary_variables, :ViscosityFactor) && !haskey(iso.secondary_variables, :ConductivityFactor)
end

@testset "Factors scale the viscosity, conductivity and diffusivity" begin
    # Factors of 2 on the viscosity, 3 on the conductivity and 1/2 on the
    # diffusivity in every cell give what those parameters give
    k = Mocca.convert_constants(haghpanah_constants())
    ads(D) = Adsorbent(; isotherm = k.adsorbent.isotherm, density = k.adsorbent.density, heat_capacity = k.adsorbent.heat_capacity,
        particle_diameter = k.adsorbent.particle_diameter, particle_porosity = k.adsorbent.particle_porosity,
        tortuosity = k.adsorbent.tortuosity, adsorbed_heat_capacity = k.adsorbent.adsorbed_heat_capacity,
        mass_transfer = LinearDrivingForce(D, k.adsorbent.tortuosity, k.adsorbent.particle_porosity, k.adsorbent.particle_diameter))
    D = k.gas.diffusivity
    function run(transport, adsorbent; parameters = Dict{Symbol, Any}())
        fs = Flowsheet(k.gas, k.feed)
        add_bed!(fs, :Bed, adsorbent, k.column; ncells = 8, dispersion = 1e-3)
        add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
        add_device!(fs, :V_product, :product; at = (:Bed, :top))
        if !isnothing(transport)
            m = fs[:Bed]
            s = m.system
            sys = FixedBed(; isotherm = s.isotherm, mass_transfer = s.mass_transfer, molecular_masses = s.molecular_masses,
                component_names = s.component_names, thermodynamics = s.thermodynamics,
                adsorbed_thermodynamics = s.adsorbed_thermodynamics, transport = transport, thermal = s.thermal)
            Mocca.network(fs).units[:Bed] = Mocca.setup_process_model(sys, m.data_domain)
        end
        p = Process(fs, Cycle([Adsorb(40.0; velocity = k.velocity)]); timestepping = TimeStepping(max_dt = 2.0, first_dt = 0.1),
            initial_state = Dict(:Bed => (pressure = k.p_high, temperature = 298.15)), parameters = parameters)
        return simulate(p).states[end][:Bed]
    end
    a = run(TransportReference(_FixedTransport(), 0.5, 1 / 3, 2.0), ads(D))
    ref = Mocca.fixed_bed_model(k.gas, k.feed, k.adsorbent, k.column; ncells = 8, dispersion = 1e-3)
    prm = Jutul.setup_parameters(ref)
    b = run(nothing, ads(D / 2); parameters = Dict(:Bed => Dict(:FluidViscosity => 2 .* prm[:FluidViscosity],
        :ThermalConductivities => 3 .* prm[:ThermalConductivities])))
    c = run(nothing, ads(D))
    for key in (:Pressure, :Temperature, :y, :AdsorbedConcentration)
        @test a[key] ≈ b[key] rtol = 1e-6
    end
    # And they matter
    @test !isapprox(a[:AdsorbedConcentration], c[:AdsorbedConcentration], rtol = 1e-3)
end

@testset "Gradients with transport properties that change" begin
    p = four_stage_vsa_process(ncells = 6, durations = (10.0, 15.0, 20.0, 25.0), P_H = 1.1bar,
        timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05), gas = Mocca.haghpanah_2013_gas(transport = :sutherland))
    add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.3bar))
    add_variable!(p, :rho, (:parameter, :Bed, :AdsorbentDensity); bounds = (1000.0, 1300.0))
    g = performance_gradient(p; progress = false)
    @test g.steady_state.converged
    value() = performance(cyclic_steady_state(p; tol = 1e-11, progress = false))
    for n in (:P_L, :rho)
        x0 = p[n]
        h = 1e-4 * x0
        p[n] = x0 + h; up = value()
        p[n] = x0 - h; down = value()
        p[n] = x0
        for m in (:purity, :recovery, :specific_energy)
            @test g.gradients[m][n] ≈ (up[m] - down[m]) / 2h rtol = 1e-4
        end
    end
end
