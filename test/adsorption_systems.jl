@testset "AdsorptionSystem Construction" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    system = haghpanah_bed(constants; ncells = 3).system

    @test system isa Mocca.AdsorptionSystem
    @test Mocca.number_of_components(system) == 2
    @test system.component_names == ["CO2", "N2"]
    @test system.isotherm isa Mocca.DualSiteLangmuir
    @test system.mass_transfer isa Mocca.LinearDrivingForce
    @test system.molecular_masses == SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2)
    # The heat capacities per mass of the constants, as molar heat capacities
    M = system.molecular_masses
    @test system.thermodynamics isa Mocca.IdealGasConstantCp
    @test Mocca.component_heat_capacities(system.thermodynamics, 350.0) ≈ constants.C_pg .* M
    @test Mocca.component_heat_capacities(system.adsorbed_thermodynamics, 350.0) ≈ constants.C_pa .* M
end

@testset "AdsorptionSystem Properties" begin
    system = haghpanah_bed(ncells = 3).system

    @test Mocca.component_names(system) == ["CO2", "N2"]
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
    @test Mocca.number_of_components(system) == 3
    @test system.component_names == ["CO2", "N2", "H2O"]
    @test length(system.molecular_masses) == 3
    # Heat capacities per mass make constant molar heat capacities, and the
    # adsorbed phase has those of the gas
    @test system.thermodynamics.cp ≈ [846.0 * 44.01e-3, 1040.0 * 28.0e-3, 1996.0 * 18.015e-3]
    @test system.adsorbed_thermodynamics == system.thermodynamics
    shomate = Mocca.IdealGasShomate([(24.99735, 55.18696, -33.69137, 7.948387, -0.136638),
        (28.98641, 1.853978, -9.647459, 16.63537, 0.000117), (30.092, 6.832514, 6.793435, -2.53448, 0.082139)])
    s2 = Mocca.FixedBed(isotherm = isotherm, mass_transfer = mass_transfer, molecular_masses = SVector(44.01e-3, 28.0e-3, 18.015e-3),
        component_names = ["CO2", "N2", "H2O"], thermodynamics = shomate)
    @test s2.thermodynamics === s2.adsorbed_thermodynamics === shomate
    @test isbits(s2.thermodynamics)
    @test_throws ErrorException Mocca.FixedBed(isotherm = isotherm, mass_transfer = mass_transfer, molecular_masses = SVector(44.01e-3, 28.0e-3, 18.015e-3),
        component_names = ["CO2", "N2", "H2O"], thermodynamics = shomate, heat_capacity_gas = SVector(846.0, 1040.0, 1996.0))
end
