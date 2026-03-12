@testset "Mass Transfer Rate" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.TwoComponentAdsorptionSystem(constants)

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

    # LinearDrivingForce is constructed in the system
    @test system.mass_transfer isa Mocca.LinearDrivingForce
end
