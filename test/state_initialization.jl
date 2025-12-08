@testset "State Initialization" begin
    # Set up system and model
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);
    # Create a simple mesh
    ncells = 10
    model = Mocca.setup_process_model(system; ncells = ncells);

    # Test state initialization
    bar = 1e5  # Pa
    P_init = 1.0 * bar
    T_init = 298.15  # K
    Tw_init = constants.T_a

    # Initial composition (very small CO2, mostly N2)
    yCO2 = fill(1e-10, ncells)
    y_init = hcat(yCO2, 1 .- yCO2)'

    state0 = Mocca.setup_process_state(model;
        Pressure = P_init,
        Temperature = T_init,
        WallTemperature = Tw_init,
        y = y_init
    )
    parameters = Mocca.setup_process_parameters(model);

    @test state0 isa Dict
    @test parameters isa Dict

    # Check that all required fields are present
    @test haskey(state0, :Pressure)
    @test haskey(state0, :y)
    @test haskey(state0, :AdsorbedConcentration)
    @test haskey(state0, :Temperature)
    @test haskey(state0, :WallTemperature)

    # Check dimensions
    @test length(state0[:Pressure]) == ncells
    @test size(state0[:y]) == (2, ncells)  # 2 components, ncells cells
    @test size(state0[:AdsorbedConcentration]) == (2, ncells)
    @test length(state0[:Temperature]) == ncells
    @test length(state0[:WallTemperature]) == ncells

    # Check initial values
    @test all(state0[:Pressure] .≈ P_init)
    @test all(state0[:Temperature] .≈ T_init)
    @test all(state0[:WallTemperature] .≈ Tw_init)

    # Check that compositions are normalized
    for i in 1:ncells
        @test sum(state0[:y][:, i]) ≈ 1.0 atol=1e-10
    end

    # Check that initial CO2 concentration is very small
    @test all(state0[:y][1, :] .≈ 1e-10)
    @test all(state0[:y][2, :] .≈ 1.0 - 1e-10)

    # Check that parameters contain expected fields
    @test haskey(parameters, :SolidVolume)
    @test haskey(parameters, :FluidVolume)
    @test length(parameters[:SolidVolume]) == ncells
    @test length(parameters[:FluidVolume]) == ncells

end

@testset "State Initialization with different size arrays" begin
    # Set up system and model
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);
    # Create a simple mesh
    ncells = 10
    model = Mocca.setup_process_model(system; ncells = ncells);

    # Test state initialization
    bar = 1e5  # Pa
    P_init = 1.0 * bar
    T_init = 298.15  # K
    Tw_init = constants.T_a

    # Initial composition (very small CO2, mostly N2)
    yCO2 = 1e-10
    y_init = [yCO2, 1 .- yCO2]

    state0 = Mocca.setup_process_state(model;
        Pressure = P_init,
        Temperature = T_init,
        WallTemperature = Tw_init,
        y = y_init
    )

 # Check that all required fields are present
    @test haskey(state0, :Pressure)
    @test haskey(state0, :y)
    @test haskey(state0, :AdsorbedConcentration)
    @test haskey(state0, :Temperature)
    @test haskey(state0, :WallTemperature)

    # Check dimensions
    @test length(state0[:Pressure]) == ncells
    @test size(state0[:y]) == (2, ncells)  # 2 components, ncells cells
    @test size(state0[:AdsorbedConcentration]) == (2, ncells)
    @test length(state0[:Temperature]) == ncells
    @test length(state0[:WallTemperature]) == ncells

end


@testset "State Initialization with Different Conditions" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    system = Mocca.TwoComponentAdsorptionSystem(constants);

    ncells = 5
    model = Mocca.setup_process_model(system; ncells = ncells);

    # Test with different pressures
    bar = 1e5
    P_high = 5.0 * bar
    P_low = 0.5 * bar
    T_init = 298.15
    Tw_init = constants.T_a

    yCO2 = fill(0.1, ncells)  # 10% CO2
    y_init = hcat(yCO2, 1 .- yCO2)'

    state_high = Mocca.setup_process_state(model;
        Pressure = P_high,
        Temperature = T_init,
        WallTemperature = Tw_init,
        y = y_init
    )
    state_low = Mocca.setup_process_state(model;
        Pressure = P_low,
        Temperature = T_init,
        WallTemperature = Tw_init,
        y = y_init
    )

    @test all(state_high[:Pressure] .≈ P_high)
    @test all(state_low[:Pressure] .≈ P_low)

    # Test with different temperatures
    T_high = 350.0
    T_low = 250.0
    P_init = 1.0 * bar

    state_hot = Mocca.setup_process_state(model;
        Pressure = P_init,
        Temperature = T_high,
        WallTemperature = Tw_init,
        y = y_init
    )
    state_cold = Mocca.setup_process_state(model;
        Pressure = P_init,
        Temperature = T_low,
        WallTemperature = Tw_init,
        y = y_init
    )

    @test all(state_hot[:Temperature] .≈ T_high)
    @test all(state_cold[:Temperature] .≈ T_low)
end
