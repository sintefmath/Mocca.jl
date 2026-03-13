using Test
using Mocca
using Jutul
using JutulDarcy
using LinearAlgebra
using StaticArrays

@testset "Mocca.jl Tests" begin
    @testset "Input Parsing" begin
        include("data_input.jl")
    end
    @testset "Adsorption Systems" begin
        include("adsorption_systems.jl")
    end

    @testset "Mass Transfer" begin
        include("mass_transfer.jl")
    end

    @testset "Permeability and Dispersion" begin
        include("permeability.jl")
    end

    @testset "State Initialization" begin
        include("state_initialization.jl")
    end

    @testset "Isotherm Models" begin
        include("isotherms.jl")
    end

    @testset "Integration with Jutul" begin
        include("jutul_integration.jl")
    end

    @testset "Regression Tests" begin
        include("regression/dcb_regression.jl")
        include("regression/cyclic_regression.jl")
    end

end
