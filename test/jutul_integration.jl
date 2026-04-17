@testset "Jutul Integration Tests" begin
    # Test that Mocca systems work with Jutul framework
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.AdsorptionSystem(constants)
    # Create a simple mesh
    ncells = 10
    model = Mocca.setup_process_model(system, constants; ncells = ncells)

    domain = model.data_domain
    mesh = domain.representation

    @test model.data_domain.representation isa Jutul.CartesianMesh
    @test Jutul.number_of_cells(mesh) == ncells

    # Test domain setup
    @test domain isa Jutul.DataDomain
    @test haskey(domain, :porosity)
    @test haskey(domain, :permeability)
    @test haskey(domain, :dx)

    # Column-entity data
    @test haskey(domain, :diffusion_coefficient, Mocca.Column())
    @test haskey(domain, :thermal_conductivity, Mocca.Column())
    @test haskey(domain, :r_in, Mocca.Column())
    @test haskey(domain, :r_out, Mocca.Column())

    # Test model creation
    model2 = Jutul.SimulationModel(domain, system)

    @test model2 isa Jutul.SimulationModel
    @test model2.system === system
    @test model2.data_domain === domain

    # Test that the model has the expected properties
    @test Jutul.number_of_cells(model2.domain) == ncells
    @test Mocca.number_of_components(model2.system) == 2
end

@testset "Model State Setup" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.AdsorptionSystem(constants)
    # Create a simple mesh
    ncells = 3
    model = Mocca.setup_process_model(system, constants; ncells = ncells)

    domain = model.data_domain

    P_init = 1e5
    T_init = 298.15
    Tw_init = constants.T_a
    yCO2 = 0.05
    y_init = [yCO2, 1 .- yCO2]

    state = Mocca.setup_process_state(model;
        Pressure = P_init,
        Temperature = T_init,
        WallTemperature = Tw_init,
        y = y_init
    )

    parameters = Mocca.setup_process_parameters(model)

    # Test that state is compatible with model
    @test Jutul.number_of_cells(domain) == length(state[:Pressure])

    # Test parameter setup
    @test parameters isa Dict
    @test haskey(parameters, :SolidVolume)
    @test haskey(parameters, :FluidVolume)

    # Column-entity parameters should be present
    @test haskey(parameters, :AdsorbentDensity)
    @test haskey(parameters, :FluidViscosity)
    @test haskey(parameters, :WallCrossSectionArea)
end

@testset "Domain Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.AdsorptionSystem(constants)
    ncells = 4
    model = Mocca.setup_process_model(system, constants; ncells = ncells)

    domain = model.data_domain
    mesh = domain.representation

    # Test domain properties are correctly set
    perm = Mocca.compute_permeability(constants.Φ, constants.d_p)
    disp = Mocca.compute_dispersion(constants.D_m, constants.V0_inter, constants.d_p)

    @test all(domain[:porosity] .== constants.Φ)
    @test all(domain[:permeability] .== perm)

    # Column-entity checks
    @test first(domain[:diffusion_coefficient, Mocca.Column()]) == disp
    @test first(domain[:thermal_conductivity, Mocca.Column()]) == constants.K_z
    @test first(domain[:r_in, Mocca.Column()]) == constants.r_in
    @test first(domain[:r_out, Mocca.Column()]) == constants.r_out

    # Test dx calculation
    expected_dx = constants.L / ncells
    @test all(domain[:dx] .≈ expected_dx)

    # Test volumes
    dr = sqrt(pi*constants.r_in^2)
    expected_volume = dr^2 * expected_dx  # Cell volume
    @test all(domain[:volumes] .≈ expected_volume)
end
