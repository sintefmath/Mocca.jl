using Test
using Mocca
using Jutul
using JutulDarcy
using MAT
using LinearAlgebra
using StaticArrays

@testset "Mocca.jl Tests" begin

    @testset "Adsorption Systems" begin
        include("adsorption_systems.jl")
    end

    @testset "Equilibrium Calculations" begin
        include("equilibrium.jl")
    end

    @testset "Permeability and Dispersion" begin
        include("permeability.jl")
    end

    @testset "State Initialization" begin
        include("state_initialization.jl")
    end

    @testset "Integration with Jutul" begin
        include("jutul_integration.jl")
    end

end
