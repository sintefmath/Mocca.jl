@testset "AdsorptionSystem Construction" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.AdsorptionSystem(constants)

    @test system isa Mocca.AdsorptionSystem
    @test JutulDarcy.number_of_components(system) == 2
    @test JutulDarcy.number_of_phases(system) == 1
    @test JutulDarcy.get_reference_phase_index(system) == 1
    @test system.component_names == ["CO2", "N2"]
    @test system.isotherm isa Mocca.DualSiteLangmuir
    @test system.mass_transfer isa Mocca.LinearDrivingForce
    @test system.molecular_masses == SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2)
    @test system.heat_capacity_gas == constants.C_pg
    @test system.heat_capacity_adsorbed == constants.C_pa
end

@testset "AdsorptionSystem Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.AdsorptionSystem(constants)

    @test !JutulDarcy.has_other_phase(system)
    @test JutulDarcy.phase_names(system) == ["gas"]
    @test collect(JutulDarcy.eachphase(system)) == [1]
end

@testset "N-component AdsorptionSystem" begin
    # Test construction of a 3-component adsorption column system
    # Note that the H2O parameters are synthetic and chosen only to
    # exercise the 3-component code path. They are not physically accurate.
    # TODO: Replace with physically accurate parameters or add an example.
    isotherm = Mocca.DualSiteLangmuir(
        qsb = SVector(3489.44, 6613.55, 1000.0),
        b0  = SVector(8.65e-7, 2.5e-6, 1.0e-6),
        ΔUb = SVector(-36641.21, -15800.0, -20000.0),
        qsd = SVector(2872.35, 0.0, 500.0),
        d0  = SVector(2.63e-8, 0.0, 1.0e-8),
        ΔUd = SVector(-35690.66, 0.0, -10000.0),
    )

    mass_transfer = Mocca.LinearDrivingForce(
        D_m = 1.6e-5, τ = 3.0, ϵ_p = 0.35, d_p = 2e-3,
    )

    system = Mocca.AdsorptionSystem(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = SVector(44.01e-3, 28.0e-3, 18.015e-3),
        component_names = ["CO2", "N2", "H2O"],
        heat_capacity_gas = SVector(846.0, 1040.0, 1996.0),
        heat_capacity_adsorbed = SVector(846.0, 1040.0, 1996.0),
    )

    @test system isa Mocca.AdsorptionSystem{3}
    @test JutulDarcy.number_of_components(system) == 3
    @test system.component_names == ["CO2", "N2", "H2O"]
    @test length(system.molecular_masses) == 3
end
