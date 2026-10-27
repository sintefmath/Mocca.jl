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
Mocca.molecular_diffusivities(::_FixedTransport, y, T, P) = SVector(one(T / P), one(T / P))

@testset "Chapman–Enskog transport properties" begin
    gas = Gas(["CO2", "N2", "H2O"]; transport = :chapman_enskog)
    t = gas.transport
    @test t isa ChapmanEnskogTransport && isbits(t)
    # NIST Chemistry WebBook at 300 and 400 K and 0.01 bar: within 3% for CO2
    # and N2; water, a polar molecule, is 10% high at 300 K
    μ(T) = Mocca.Physics.component_viscosities(t, T)
    @test μ(300.0)[1:2] ≈ [14.994e-6, 17.877e-6] rtol = 0.03
    @test μ(400.0)[1:2] ≈ [19.621e-6, 22.200e-6] rtol = 0.03
    @test μ(300.0)[3] ≈ 9.7659e-6 rtol = 0.11
    @test μ(400.0)[3] ≈ 13.354e-6 rtol = 0.07
    # A pure component has its own viscosity, the mixture lies between
    @test mixture_viscosity(t, [1.0, 0.0, 0.0], 350.0) ≈ μ(350.0)[1]
    @test μ(300.0)[3] < mixture_viscosity(t, [0.15, 0.82, 0.03], 300.0) < μ(300.0)[2]

    # CO2 in N2 at 25 °C and 1 atm: about 1.6e-5 m²/s measured
    D12(T, P) = Mocca.Physics.binary_diffusivity(t, 1, 2, T, P)
    @test D12(298.15, 101325.0) ≈ 1.65e-5 rtol = 0.07
    @test D12(298.15, 101325.0) ≈ Mocca.Physics.binary_diffusivity(t, 2, 1, 298.15, 101325.0)
    @test D12(350.0, 0.5e5) / D12(350.0, 1e5) ≈ 2
    # In a binary mixture each component has the binary diffusivity, whatever
    # the composition
    b = Gas(["CO2", "N2"]; transport = :chapman_enskog).transport
    for y in ([0.15, 0.85], [0.7, 0.3], [1 - 1e-10, 1e-10])
        @test molecular_diffusivities(b, y, 320.0, 1e5) ≈ fill(Mocca.Physics.binary_diffusivity(b, 1, 2, 320.0, 1e5), 2) rtol = 1e-9
    end
    # In nearly pure N2 with traces of CO2 and water, the diffusivity of N2 is
    # the mean of its binary ones, whatever the ratio of the traces, and its
    # derivatives stay moderate, so that Newton's method converges
    Dn(y) = molecular_diffusivities(t, y, 320.0, 1e5)[2]
    bin = [Mocca.Physics.binary_diffusivity(t, 2, j, 320.0, 1e5) for j in (1, 3)]
    for y in ([1e-10, 1 - 2e-10, 1e-10], [1e-6 * 1e-4, 1 - 1e-6, 1e-6], [0.0, 1.0, 0.0])
        @test minimum(bin) <= Dn(y) <= maximum(bin)
    end
    @test Dn([1e-12, 1 - 2e-10, 2e-10 - 1e-12]) ≈ Dn([2e-10 - 1e-12, 1 - 2e-10, 1e-12]) rtol = 1e-3
    # (of order D/1e-6, the floor on the traces, not D/1e-10)
    g = ForwardDiff.gradient(Dn, [1e-10, 1 - 2e-10, 1e-10])
    @test maximum(abs, g) < Dn([1e-10, 1 - 2e-10, 1e-10]) / 1e-6
    @test_throws ArgumentError ChapmanEnskogTransport(; sigma = [3.9, 3.8], epsilon = [195.0], molar_masses = [44e-3, 28e-3])

    # Automatic differentiation without allocating
    ys = SVector(0.15, 0.82, 0.03)
    T = ForwardDiff.Dual(330.0, 1.0)
    f(t, T) = (mixture_viscosity(t, ys, T), molecular_diffusivities(t, ys, T, 1e5))
    f(t, T)
    @test (@allocated f(t, T)) == 0
    h = 1e-3
    @test ForwardDiff.partials(f(t, T)[1])[1] ≈ (mixture_viscosity(t, ys, 330.0 + h) - mixture_viscosity(t, ys, 330.0 - h)) / 2h rtol = 1e-6
end

@testset "Factors against the reference" begin
    t = Gas(["CO2", "N2"]; transport = :chapman_enskog).transport
    y, T, P = [0.15, 0.85], 298.15, 1e5
    r = TransportReference(t, y, T, P)
    @test Mocca.viscosity_factor(r, y, T) == 1 && all(==(1), Mocca.diffusivity_factors(r, y, T, P))
    @test r.pressure == P
    @test Mocca.viscosity_factor(r, y, 350.0) ≈ mixture_viscosity(t, y, 350.0) / mixture_viscosity(t, y, T)
    @test Mocca.diffusivity_factors(r, y, T, P / 2) ≈ [2, 2]
    c = TransportReference(ConstantTransport())
    @test Mocca.is_constant_transport(c) && !Mocca.is_constant_transport(r)
    @test TransportReference(ConstantTransport(), y, T, P) == c
    @test Mocca.viscosity_factor(c, y, 400.0) === 1.0

    # The linear driving force scales the molecular diffusivity of each component
    mt = LinearDrivingForce(D_m = 1.6e-5, τ = 3.0, ϵ_p = 0.35, d_p = 2e-3)
    C, q, qs = SVector(5.0, 30.0), SVector(100.0, 50.0), SVector(800.0, 60.0)
    @test compute_mass_transfer_rate(mt, C, q, qs, SVector(2.0, 0.5)) ≈
        compute_mass_transfer_rate(LinearDrivingForce(D_m = SVector(3.2e-5, 0.8e-5), τ = 3.0, ϵ_p = 0.35, d_p = 2e-3), C, q, qs)

    # Wakao and Funazkri: ε D_ax = 20 D_m + 0.5 |u| d_p
    wf = WakaoFunazkri(; particle_diameter = 2e-3, porosity = 0.37, diffusivities = [1.6e-5, 1.6e-5])
    @test Mocca.dispersion_coefficient(wf, 1.6e-5, -0.5) ≈ (20 * 1.6e-5 + 0.5 * 0.5 * 2e-3) / 0.37
end

@testset "Gas transport" begin
    @test Gas(["CO2", "N2"]).transport === ConstantTransport()
    @test Gas(["CO2", "N2"]; transport = _FixedTransport()).transport === _FixedTransport()
    @test_throws MoccaInputError Gas(["CO2", "N2"]; transport = :sutherland)
    # No Lennard-Jones parameters for CO; a component can be given them
    @test_throws MoccaInputError Gas(["CO2", "CO"]; transport = :chapman_enskog)
    register_component!(Component("Ne"; molar_mass = 20.18e-3, heat_capacity = 1030.0, viscosity = 3.2e-5, lennard_jones = (2.820, 32.8)))
    @test Gas(["CO2", "Ne"]; transport = :chapman_enskog).transport.σ[2] == 2.820
    @test_throws MoccaInputError Component("Xx"; molar_mass = 1.0, heat_capacity = 1.0, viscosity = 1.0, lennard_jones = (2.8,))

    # The viscosity and diffusivities at the feed: given, from the model, or
    # the constant defaults
    g = Gas(["CO2", "N2"]; transport = :chapman_enskog)
    y = [0.15, 0.85]
    @test Mocca.gas_viscosity(g, y, 350.0) ≈ mixture_viscosity(g.transport, y, 350.0)
    @test Mocca.gas_viscosity(Gas(["CO2", "N2"]), y, 350.0) == Mocca.gas_viscosity(Gas(["CO2", "N2"]), y)
    @test Mocca.gas_viscosity(Gas(["CO2", "N2"]; transport = :chapman_enskog, viscosity = 1.7e-5), y, 350.0) == 1.7e-5
    @test Mocca.gas_diffusivities(g, y, 300.0, 1e5) ≈ molecular_diffusivities(g.transport, y, 300.0, 1e5)
    @test Mocca.gas_diffusivities(Gas(["CO2", "N2"]), y, 300.0, 1e5) == [1.6e-5, 1.6e-5]
    @test Mocca.gas_diffusivities(Gas(["CO2", "N2"]; diffusivity = [1e-5, 2e-5], transport = :chapman_enskog), y, 300.0, 1e5) == [1e-5, 2e-5]

    # In a process file
    p = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = :chapman_enskog))
    d = Mocca.process_to_dict(p)
    @test d["gas"]["transport"] == "chapman_enskog"
    q = read_process(write_process(joinpath(mktempdir(), "t.json"), p))
    @test q.flowsheet.gas.transport == p.flowsheet.gas.transport
    @test q.flowsheet[:Bed].system.transport == p.flowsheet[:Bed].system.transport
    @test q.flowsheet[:Bed].system.dispersion == p.flowsheet[:Bed].system.dispersion
    @test !haskey(Mocca.process_to_dict(four_stage_vsa_process(ncells = 4))["gas"], "transport")
    r = four_stage_vsa_process(ncells = 4, gas = Mocca.haghpanah_2013_gas(transport = _FixedTransport()))
    @test_throws MoccaInputError Mocca.process_to_dict(r)
    # A component of the user's own, with its data
    ne = Gas(["N2", "Ne"]; transport = :chapman_enskog)
    pn = four_stage_vsa_process(ncells = 4, gas = ne, adsorbent = Adsorbent(; isotherm = load_adsorbent("zeolite_13x_haghpanah_2013").isotherm,
        density = 1130.0, heat_capacity = 1070.0, particle_diameter = 2e-3, particle_porosity = 0.35))
    dn = Mocca.process_to_dict(pn)
    @test only(dn["components"])["lennard_jones"] == [2.820, 32.8]
end

@testset "A bed with transport properties that change" begin
    # The bed refers the viscosity and diffusivities to the feed, and has the
    # variables only when its transport changes
    gas = Mocca.haghpanah_2013_gas(transport = :chapman_enskog)
    p = four_stage_vsa_process(ncells = 4, gas = gas)
    bed = p.flowsheet[:Bed]
    feed = p.flowsheet.feed
    @test bed.system.transport == TransportReference(gas.transport, feed.composition, feed.temperature, feed.pressure)
    @test haskey(bed.secondary_variables, :ViscosityFactor) && haskey(bed.secondary_variables, :MolecularDiffusivity)
    @test bed.secondary_variables[:AdsorptionMassTransfer] isa Mocca.ScaledAdsorptionMassTransfer
    # Dispersion by Wakao and Funazkri from the diffusivities at the feed (the
    # preset gives 1.6e-5 m²/s), scaled by a parameter
    @test bed.system.dispersion == WakaoFunazkri(; particle_diameter = 2e-3, porosity = 0.37, diffusivities = [1.6e-5, 1.6e-5])
    @test only(p.parameters[:Bed][:DispersionScale]) == 1.0
    plain = four_stage_vsa_process(ncells = 4).flowsheet[:Bed]
    @test !haskey(plain.secondary_variables, :ViscosityFactor) && !haskey(plain.secondary_variables, :MolecularDiffusivity)
    @test plain.secondary_variables[:AdsorptionMassTransfer] isa Mocca.AdsorptionMassTransfer
    @test plain.system.dispersion === ConstantDispersion()
    # A dispersion given stays constant
    fixed = Mocca.fixed_bed_model(gas, feed, load_adsorbent("zeolite_13x_haghpanah_2013"), Mocca.haghpanah_2013_column();
        ncells = 3, dispersion = 1e-3)
    @test fixed.system.dispersion === ConstantDispersion() && !haskey(fixed.secondary_variables, :MolecularDiffusivity)
end

# A bed of the Mocca 0.1.0 constants with its transport and dispersion
# replaced, simulated through a 40 s breakthrough
function _transport_breakthrough(transport, dispersion, adsorbent, k; parameters = Dict{Symbol, Any}(), D_L = 1e-3)
    fs = Flowsheet(k.gas, k.feed)
    add_bed!(fs, :Bed, adsorbent, k.column; ncells = 8, dispersion = D_L)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    if !isnothing(transport)
        m = fs[:Bed]
        s = m.system
        sys = FixedBed(; isotherm = s.isotherm, mass_transfer = s.mass_transfer, molecular_masses = s.molecular_masses,
            component_names = s.component_names, thermodynamics = s.thermodynamics,
            adsorbed_thermodynamics = s.adsorbed_thermodynamics, transport = transport, dispersion = dispersion, thermal = s.thermal)
        Mocca.network(fs).units[:Bed] = Mocca.setup_process_model(sys, m.data_domain)
    end
    p = Process(fs, Cycle([Adsorb(40.0; velocity = k.velocity)]); timestepping = TimeStepping(max_dt = 2.0, first_dt = 0.1),
        initial_state = Dict(:Bed => (pressure = k.p_high, temperature = 298.15)), parameters = parameters)
    return simulate(p).states[end][:Bed]
end

@testset "Factors scale the viscosity, diffusivity and dispersion" begin
    k = Mocca.convert_constants(haghpanah_constants())
    ads(D) = Adsorbent(; isotherm = k.adsorbent.isotherm, density = k.adsorbent.density, heat_capacity = k.adsorbent.heat_capacity,
        particle_diameter = k.adsorbent.particle_diameter, particle_porosity = k.adsorbent.particle_porosity,
        tortuosity = k.adsorbent.tortuosity, adsorbed_heat_capacity = k.adsorbent.adsorbed_heat_capacity,
        mass_transfer = LinearDrivingForce(D, k.adsorbent.tortuosity, k.adsorbent.particle_porosity, k.adsorbent.particle_diameter))
    D = k.gas.diffusivity
    compare(a, b) = for key in (:Pressure, :Temperature, :y, :AdsorbedConcentration)
        @test a[key] ≈ b[key] rtol = 1e-6
    end
    ref = Mocca.fixed_bed_model(k.gas, k.feed, k.adsorbent, k.column; ncells = 8, dispersion = 1e-3)
    prm = Jutul.setup_parameters(ref)
    # Factors of 2 on the viscosity and 1/2 on the diffusivity in every cell
    # give what those values give
    a = _transport_breakthrough(TransportReference(_FixedTransport(), 0.5, SVector(2.0, 2.0), 1e5), ConstantDispersion(), ads(D), k)
    b = _transport_breakthrough(nothing, nothing, ads(D / 2), k; parameters = Dict(:Bed => Dict(:FluidViscosity => 2 .* prm[:FluidViscosity])))
    compare(a, b)
    # And they matter
    c = _transport_breakthrough(nothing, nothing, ads(D), k)
    @test !isapprox(a[:AdsorbedConcentration], c[:AdsorbedConcentration], rtol = 1e-3)

    # Wakao and Funazkri without the velocity term (d_p = 0) and with fixed
    # diffusivities is a constant dispersion of 20 D/ε, times its scale
    unit = TransportReference(_FixedTransport(), 1.0, SVector(1.0, 1.0), 1e5)
    ε = k.column.porosity
    wf = WakaoFunazkri(; particle_diameter = 0.0, porosity = ε, diffusivities = [D, D])
    compare(_transport_breakthrough(unit, wf, ads(D), k), _transport_breakthrough(nothing, nothing, ads(D), k; D_L = 20D / ε))
    compare(_transport_breakthrough(unit, wf, ads(D), k; parameters = Dict(:Bed => Dict(:DispersionScale => [0.5]))),
        _transport_breakthrough(nothing, nothing, ads(D), k; D_L = 10D / ε))
end

@testset "Gradients with transport properties that change" begin
    p = four_stage_vsa_process(ncells = 6, durations = (10.0, 15.0, 20.0, 25.0), P_H = 1.1bar,
        timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05), gas = Mocca.haghpanah_2013_gas(transport = :chapman_enskog))
    add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.3bar))
    add_variable!(p, :rho, (:parameter, :Bed, :AdsorbentDensity); bounds = (1000.0, 1300.0))
    add_variable!(p, :disp, (:parameter, :Bed, :DispersionScale); bounds = (0.5, 2.0))
    g = performance_gradient(p; progress = false)
    @test g.steady_state.converged
    value() = performance(cyclic_steady_state(p; tol = 1e-11, progress = false))
    for n in (:P_L, :rho, :disp)
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
