@testset "Jutul Integration Tests" begin
    # Test that Mocca systems work with Jutul framework
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);
    # Create a simple mesh
    ncells = 10
    model = Mocca.setup_process_model(system; ncells = ncells);

    domain = model.data_domain
    mesh = domain.representation

    @test model.data_domain.representation isa Jutul.CartesianMesh
    @test Jutul.number_of_cells(mesh) == ncells

    # Test domain setup
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
  # Test that Mocca systems work with Jutul framework
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);
    # Create a simple mesh
    ncells = 3
    model = Mocca.setup_process_model(system; ncells = ncells);

    domain = model.data_domain;
    mesh = domain.representation;

    # Test state setup using Jutul's setup_state
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

    parameters = Mocca.setup_process_parameters(model);
    
    # Test that state is compatible with model
    @test Jutul.number_of_cells(domain) == length(state[:Pressure])

    # Test parameter setup
    @test parameters isa Dict
    @test haskey(parameters, :SolidVolume)
    @test haskey(parameters, :FluidVolume)
end

@testset "Domain Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);

    ncells = 4
    model = Mocca.setup_process_model(system; ncells = ncells);

    domain = model.data_domain
    mesh = domain.representation

    # Test domain properties are correctly set
    @test all(domain[:porosity] .== system.p.Φ)
    @test all(domain[:permeability] .== system.permeability)
    @test all(domain[:diffusion_coefficient] .== system.dispersion)
    @test all(domain[:thermal_conductivity] .== system.p.K_z)

    # Test dx calculation
    expected_dx = constants.L / ncells
    @test all(domain[:dx] .≈ expected_dx)

    # Test volumes
    dr = sqrt(pi*constants.r_in^2)
    expected_volume =dr.^2 * expected_dx  # Cell volume
    @test all(domain[:volumes] .≈ expected_volume)
end
