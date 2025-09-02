@testset "AdsorptionSystem Construction" begin
    # Test TwoComponentAdsorptionSystem creation
    constants = Mocca.HaghpanahConstants{Float64}()
    permeability = Mocca.compute_permeability(constants)
    dispersion = Mocca.calc_dispersion(constants)

    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = permeability,
        dispersion = dispersion,
        p = constants
    )

    @test system isa Mocca.TwoComponentAdsorptionSystem
    @test JutulDarcy.number_of_components(system) == 2
    @test JutulDarcy.number_of_phases(system) == 1
    @test JutulDarcy.get_reference_phase_index(system) == 1
    @test system.component_names == ["CO2", "N2"]
    @test system.permeability == permeability
    @test system.dispersion == dispersion
    @test system.p === constants
end

@testset "AdsorptionSystem Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    permeability = Mocca.compute_permeability(constants)
    dispersion = Mocca.calc_dispersion(constants)

    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = permeability,
        dispersion = dispersion,
        p = constants
    )

    @test !JutulDarcy.has_other_phase(system)
    @test JutulDarcy.phase_names(system) == ["gas"]
    @test collect(JutulDarcy.eachphase(system)) == [1]
end
