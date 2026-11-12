@testset "Mass Transfer Rate" begin
    system = haghpanah_bed(ncells = 3).system

    temperature = 298.15
    concentration = [0.1, 0.9]
    qstar = Mocca.compute_equilibrium(system.isotherm, concentration, temperature)

    # Positive rate when q < q*
    q_zero = [0.0, 0.0]
    rate = Mocca.compute_mass_transfer_rate(system.mass_transfer, concentration, q_zero, qstar)
    @test length(rate) == 2
    @test all(rate .> 0)

    # Zero rate at equilibrium
    rate_eq = Mocca.compute_mass_transfer_rate(system.mass_transfer, concentration, qstar, qstar)
    @test all(rate_eq .≈ 0.0)

    # Numerical stability with near-zero qstar
    qstar_small = [1e-10, 1e-10]
    rate_small = Mocca.compute_mass_transfer_rate(system.mass_transfer, concentration, q_zero, qstar_small)
    @test all(isfinite.(rate_small))

    # A component the sorbent does not take up (q* = 0) has no mass transfer
    rate_none = Mocca.compute_mass_transfer_rate(system.mass_transfer, concentration, q_zero, [qstar[1], 0.0])
    @test all(isfinite, rate_none) && rate_none[1] > 0 && rate_none[2] == 0

    # LinearDrivingForce is constructed in the system
    @test system.mass_transfer isa Mocca.LinearDrivingForce
end

@testset "Knudsen diffusion" begin
    M = [44.01e-3, 28.0134e-3]
    plain = Mocca.LinearDrivingForce(1.6e-5, 3.0, 0.5, 7e-4)
    kn = Mocca.LinearDrivingForce(1.6e-5, 3.0, 0.5, 7e-4; macropore_radius = 5.25e-7, molar_masses = M)
    @test !Mocca.has_knudsen(plain) && Mocca.has_knudsen(kn)
    C, q, qstar = [5.0, 35.0], [0.0, 0.0], [300.0, 200.0]
    T = 303.15
    # Without macropores the temperature does not matter
    @test Mocca.compute_mass_transfer_rate(plain, C, q, qstar, 400.0) == Mocca.compute_mass_transfer_rate(plain, C, q, qstar)
    # Molecular and Knudsen diffusion in series
    D_K = 2 / 3 .* 5.25e-7 .* sqrt.(8 * Mocca.GAS_CONSTANT * T ./ (π .* M))
    D = 1 ./ (1 / 1.6e-5 .+ 1 ./ D_K)
    @test Mocca.pore_diffusivity(kn, T) ≈ D
    @test Mocca.compute_mass_transfer_rate(kn, C, q, qstar, T) ≈ Mocca.compute_mass_transfer_rate(plain, C, q, qstar) .* D ./ 1.6e-5
    # Knudsen diffusion grows with temperature, and is faster for the lighter N2
    @test all(Mocca.pore_diffusivity(kn, 450.0) .> Mocca.pore_diffusivity(kn, T))
    @test Mocca.pore_diffusivity(kn, T)[2] > Mocca.pore_diffusivity(kn, T)[1]
    # Wide pores leave molecular diffusion
    wide = Mocca.LinearDrivingForce(1.6e-5, 3.0, 0.5, 7e-4; macropore_radius = 1e-2, molar_masses = M)
    @test Mocca.pore_diffusivity(wide, T) ≈ [1.6e-5, 1.6e-5] rtol = 1e-4
    # One molecular diffusivity per component
    per = Mocca.LinearDrivingForce([1.6e-5, 2e-5], 3.0, 0.5, 7e-4; macropore_radius = 5.25e-7, molar_masses = M)
    @test Mocca.pore_diffusivity(per, T)[2] ≈ 1 / (1 / 2e-5 + 1 / D_K[2])
    # The keyword form
    kw = Mocca.LinearDrivingForce(D_m = 1.6e-5, τ = 3.0, ϵ_p = 0.5, d_p = 7e-4, macropore_radius = 5.25e-7, molar_masses = M)
    @test Mocca.pore_diffusivity(kw, T) == Mocca.pore_diffusivity(kn, T)
    # The derivative in temperature
    f(x) = Mocca.compute_mass_transfer_rate(kn, C, q, qstar, x)[1]
    @test Jutul.ForwardDiff.derivative(f, T) ≈ (f(T + 1e-3) - f(T - 1e-3)) / 2e-3 rtol = 1e-6
    # Knudsen diffusion needs the temperature and the molar masses
    @test_throws ErrorException Mocca.compute_mass_transfer_rate(kn, C, q, qstar)
    @test_throws ErrorException Mocca.LinearDrivingForce(1.6e-5, 3.0, 0.5, 7e-4; macropore_radius = 5e-7)
    @test_throws ErrorException Mocca.LinearDrivingForce([1.6e-5, 2e-5], 3.0, 0.5, 7e-4; macropore_radius = 5e-7, molar_masses = [0.044])
    @test_throws ErrorException Mocca.LinearDrivingForce(1.6e-5, 3.0, 0.5, 7e-4; macropore_radius = -1.0)
end
