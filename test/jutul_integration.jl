@testset "Jutul Integration Tests" begin
    # Test that Mocca systems work with Jutul framework
    constants = Mocca.HaghpanahConstants{Float64}()
    permeability = Mocca.compute_permeability(constants)
    dispersion = Mocca.calc_dispersion(constants)

    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = permeability,
        dispersion = dispersion,
        p = constants
    )

    # Test domain creation
    ncells = 5
    dx = sqrt(pi * constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))

    @test mesh isa Jutul.CartesianMesh
    @test Jutul.number_of_cells(mesh) == ncells

    # Test domain setup
    domain = Mocca.mocca_domain(mesh, system)

    @test domain isa Jutul.DataDomain
    @test haskey(domain, :porosity)
    @test haskey(domain, :permeability)
    @test haskey(domain, :diffusion_coefficient)
    @test haskey(domain, :thermal_conductivity)
    @test haskey(domain, :dx)

    # Test model creation
    model = Jutul.SimulationModel(domain, system)

    @test model isa Jutul.SimulationModel
    @test model.system === system
    @test model.data_domain === domain

    # Test that the model has the expected properties
    @test Jutul.number_of_cells(model.domain) == ncells
    @test JutulDarcy.number_of_components(model.system) == 2
    @test JutulDarcy.number_of_phases(model.system) == 1
end

@testset "Model State Setup" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    permeability = Mocca.compute_permeability(constants)
    dispersion = Mocca.calc_dispersion(constants)

    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = permeability,
        dispersion = dispersion,
        p = constants
    )

    ncells = 3
    dx = sqrt(pi * constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = Mocca.mocca_domain(mesh, system)
    model = Jutul.SimulationModel(domain, system)

    # Test state setup using Jutul's setup_state
    P_init = 1e5
    T_init = 298.15
    Tw_init = constants.T_a
    yCO2 = fill(0.05, ncells)
    y_init = hcat(yCO2, 1 .- yCO2)

    state, parameters = Mocca.initialise_state_AdsorptionColumn(
        P_init, T_init, Tw_init, y_init, model
    )

    # Test that state is compatible with model
    @test Jutul.number_of_cells(model.domain) == length(state[:Pressure])

    # Test parameter setup
    @test parameters isa Dict
    @test haskey(parameters, :SolidVolume)
    @test haskey(parameters, :FluidVolume)
end

@testset "Domain Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    permeability = Mocca.compute_permeability(constants)
    dispersion = Mocca.calc_dispersion(constants)

    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = permeability,
        dispersion = dispersion,
        p = constants
    )

    ncells = 4
    dx = sqrt(pi * constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = Mocca.mocca_domain(mesh, system)

    # Test domain properties are correctly set
    @test all(domain[:porosity] .== system.p.Φ)
    @test all(domain[:permeability] .== system.permeability)
    @test all(domain[:diffusion_coefficient] .== system.dispersion)
    @test all(domain[:thermal_conductivity] .== system.p.K_z)

    # Test dx calculation
    expected_dx = constants.L / ncells
    @test all(domain[:dx] .≈ expected_dx)

    # Test volumes
    expected_volume = dx^2 * expected_dx  # Cell volume
    @test all(domain[:volumes] .≈ expected_volume)
end
