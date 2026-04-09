@testset "Permeability Calculation" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    permeability = Mocca.compute_permeability(constants.Φ, constants.d_p)

    @test permeability > 0
    @test permeability isa Real

    # Higher porosity should give higher permeability
    perm_high = Mocca.compute_permeability(0.5, constants.d_p)
    perm_low = Mocca.compute_permeability(0.2, constants.d_p)
    @test perm_high > perm_low

    # Larger particles should give higher permeability
    perm_large = Mocca.compute_permeability(constants.Φ, 0.002)
    perm_small = Mocca.compute_permeability(constants.Φ, 0.001)
    @test perm_large > perm_small
end

@testset "Dispersion Calculation" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    dispersion = Mocca.compute_dispersion(constants.D_m, constants.V0_inter, constants.d_p)

    @test dispersion > 0
    @test dispersion isa Real

    # Higher molecular diffusivity should give higher dispersion
    disp_high = Mocca.compute_dispersion(2e-5, constants.V0_inter, constants.d_p)
    disp_low = Mocca.compute_dispersion(1e-5, constants.V0_inter, constants.d_p)
    @test disp_high > disp_low

    # Higher velocity should give higher dispersion
    disp_high_vel = Mocca.compute_dispersion(constants.D_m, 0.2, constants.d_p)
    disp_low_vel = Mocca.compute_dispersion(constants.D_m, 0.1, constants.d_p)
    @test disp_high_vel > disp_low_vel

    # Larger particles should give higher dispersion
    disp_large = Mocca.compute_dispersion(constants.D_m, constants.V0_inter, 0.002)
    disp_small = Mocca.compute_dispersion(constants.D_m, constants.V0_inter, 0.001)
    @test disp_large > disp_small
end

@testset "Constants Structure" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    @test constants isa Mocca.ConstantsStruct
    @test constants.R > 0  # Gas constant should be positive
    @test 0 < constants.Φ < 1  # Porosity should be between 0 and 1
    @test constants.d_p > 0  # Particle diameter should be positive
    @test constants.D_m > 0  # Molecular diffusivity should be positive
    @test constants.V0_inter >= 0  # Velocity should be non-negative

    # Test custom constants
    custom_constants = Mocca.HaghpanahConstants{Float64}(
        Φ = 0.5,
        d_p = 0.002,
        D_m = 2e-5
    )

    @test custom_constants.Φ == 0.5
    @test custom_constants.d_p == 0.002
    @test custom_constants.D_m == 2e-5
end
