@testset "AdsorptionSystem Construction" begin
    # Test TwoComponentAdsorptionSystem creation
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants)

    @test system isa Mocca.TwoComponentAdsorptionSystem
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
    system = Mocca.TwoComponentAdsorptionSystem(constants)

    @test !JutulDarcy.has_other_phase(system)
    @test JutulDarcy.phase_names(system) == ["gas"]
    @test collect(JutulDarcy.eachphase(system)) == [1]
end
