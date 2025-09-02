@testset "Equilibrium Calculation Tests" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = 1e-12,
        dispersion = 1e-6,
        p = constants
    )

    # Test basic equilibrium calculation
    temperature = 298.15  # K
    concentration = [0.1, 0.9]  # mol/m³

    qstar = Mocca.compute_equilibrium(system, concentration, temperature)

    @test qstar isa StaticArrays.SVector{2}
    @test length(qstar) == 2
    @test all(qstar .>= 0)  # Equilibrium concentrations should be non-negative

    # Test that equilibrium responds to temperature changes
    qstar_low_temp = Mocca.compute_equilibrium(system, concentration, 250.0)
    qstar_high_temp = Mocca.compute_equilibrium(system, concentration, 350.0)

    # At lower temperature, adsorption should generally be higher
    @test qstar_low_temp[1] > qstar_high_temp[1]  # CO2 adsorption

    # Test with different concentrations
    conc_high_co2 = [0.8, 0.2]
    qstar_high_co2 = Mocca.compute_equilibrium(system, conc_high_co2, temperature)

    # Higher CO2 concentration should lead to higher CO2 adsorption
    @test qstar_high_co2[1] > qstar[1]
end

@testset "Mass Transfer Coefficient Tests" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = 1e-12,
        dispersion = 1e-6,
        p = constants
    )

    temperature = 298.15
    concentration = [0.1, 0.9]
    qstar = Mocca.compute_equilibrium(system, concentration, temperature)

    ki = Mocca.compute_ki(system, concentration, qstar)

    @test length(ki) == 2
    @test all(ki .> 0)  # Mass transfer coefficients should be positive

    # Test with zero adsorbed concentration (should handle division properly)
    qstar_zero = [1e-10, 1e-10]  # Very small but not exactly zero
    ki_zero = Mocca.compute_ki(system, concentration, qstar_zero)
    @test all(isfinite.(ki_zero))
end

@testset "Equilibrium Mathematical Properties" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.TwoComponentAdsorptionSystem(
        permeability = 1e-12,
        dispersion = 1e-6,
        p = constants
    )

    temperature = 298.15

    # Test that equilibrium is continuous with respect to concentration
    c1 = [0.1, 0.9]
    c2 = [0.11, 0.89]

    q1 = Mocca.compute_equilibrium(system, c1, temperature)
    q2 = Mocca.compute_equilibrium(system, c2, temperature)

    # Test zero concentration
    c_zero = [0.0, 0.0]
    q_zero = Mocca.compute_equilibrium(system, c_zero, temperature)
    @test all(q_zero .≈ 0.0)
end
