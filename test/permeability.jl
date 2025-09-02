@testset "Permeability Calculation" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    permeability = Mocca.compute_permeability(constants)

    @test permeability > 0
    @test permeability isa Real

    # Test that permeability depends on porosity
    constants_high_porosity = Mocca.HaghpanahConstants{Float64}(Φ = 0.5)
    constants_low_porosity = Mocca.HaghpanahConstants{Float64}(Φ = 0.2)

    perm_high = Mocca.compute_permeability(constants_high_porosity)
    perm_low = Mocca.compute_permeability(constants_low_porosity)

    # Higher porosity should give higher permeability
    @test perm_high > perm_low

    # Test that permeability depends on particle diameter
    constants_large_particles = Mocca.HaghpanahConstants{Float64}(d_p = 0.002)
    constants_small_particles = Mocca.HaghpanahConstants{Float64}(d_p = 0.001)

    perm_large = Mocca.compute_permeability(constants_large_particles)
    perm_small = Mocca.compute_permeability(constants_small_particles)

    # Larger particles should give higher permeability
    @test perm_large > perm_small
end

@testset "Dispersion Calculation" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    dispersion = Mocca.calc_dispersion(constants)

    @test dispersion > 0
    @test dispersion isa Real

    # Test that dispersion depends on molecular diffusivity
    constants_high_diff = Mocca.HaghpanahConstants{Float64}(D_m = 2e-5)
    constants_low_diff = Mocca.HaghpanahConstants{Float64}(D_m = 1e-5)

    disp_high = Mocca.calc_dispersion(constants_high_diff)
    disp_low = Mocca.calc_dispersion(constants_low_diff)

    # Higher molecular diffusivity should give higher dispersion
    @test disp_high > disp_low

    # Test that dispersion depends on velocity
    constants_high_vel = Mocca.HaghpanahConstants{Float64}(V0_inter = 0.2)
    constants_low_vel = Mocca.HaghpanahConstants{Float64}(V0_inter = 0.1)

    disp_high_vel = Mocca.calc_dispersion(constants_high_vel)
    disp_low_vel = Mocca.calc_dispersion(constants_low_vel)

    # Higher velocity should give higher dispersion
    @test disp_high_vel > disp_low_vel

    # Test that dispersion depends on particle diameter
    constants_large_particles = Mocca.HaghpanahConstants{Float64}(d_p = 0.002)
    constants_small_particles = Mocca.HaghpanahConstants{Float64}(d_p = 0.001)

    disp_large = Mocca.calc_dispersion(constants_large_particles)
    disp_small = Mocca.calc_dispersion(constants_small_particles)

    # Larger particles should give higher dispersion
    @test disp_large > disp_small
end

@testset "Constants Structure" begin
    # Test default constants
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
