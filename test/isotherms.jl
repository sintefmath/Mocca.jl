using Test
using Mocca
using StaticArrays

@testset "Isotherm Tests" begin
    constants = HaghpanahConstants{Float64}()
    isotherm = DualSiteLangmuir(constants)

    @testset "compute_equilibrium" begin
        T = 298.15
        C = @SVector [100.0, 200.0]

        qstar = compute_equilibrium(isotherm, C, T)

        @test qstar isa SVector{2, Float64}
        @test all(qstar .>= 0.0)  # Loadings should be non-negative
        @test qstar[1] > 0.0  # CO2 should have non-zero loading

        # Higher concentration → higher loading
        C2 = @SVector [200.0, 200.0]
        qstar2 = compute_equilibrium(isotherm, C2, T)
        @test qstar2[1] > qstar[1]

        # Higher temperature → lower loading
        @test compute_equilibrium(isotherm, C, 350.0)[1] < qstar[1]
    end

    @testset "compute_enthalpy" begin
        # Enthalpy should be independent of concentration and temperature for this model
        C = @SVector [100.0, 200.0]
        ΔH = compute_enthalpy(isotherm, C, 298.15)
        @test ΔH isa SVector{2, Float64}
        @test all(ΔH .< 0.0)  # Adsorption is exothermic

        # State-independent: different C, T give the same result
        @test compute_enthalpy(isotherm, @SVector([200.0, 100.0]), 350.0) === ΔH
    end

    @testset "Edge cases" begin
        T = 298.15
        C = @SVector [100.0, 200.0]

        # Zero concentration → zero loading
        C_zero = @SVector [0.0, 0.0]
        @test all(compute_equilibrium(isotherm, C_zero, T) .≈ 0.0)

        # High temperature → reduced but non-negative loading
        qstar_high = compute_equilibrium(isotherm, C, 1000.0)
        @test all(qstar_high .>= 0.0)
        @test qstar_high[1] < compute_equilibrium(isotherm, C, T)[1]

        # Low temperature → increased loading
        @test compute_equilibrium(isotherm, C, 250.0)[1] > compute_equilibrium(isotherm, C, T)[1]
    end
end
