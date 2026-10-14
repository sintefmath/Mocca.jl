using Test
using Mocca
using Jutul
using LinearAlgebra
using StaticArrays

include("haghpanah_cases.jl")

@testset "Mocca.jl Tests" begin
    @testset "Input Parsing" begin
        include("data_input.jl")
    end
    @testset "Adsorption Systems" begin
        include("adsorption_systems.jl")
    end

    @testset "Architecture" begin
        include("architecture.jl")
    end

    @testset "Module names" begin
        include("module_names.jl")
    end

    @testset "Flowsheets" begin
        include("flowsheet.jl")
    end

    @testset "Process metrics and cyclic steady state" begin
        include("process.jl")
    end

    @testset "Force and cyclic steady state gradients" begin
        include("gradients.jl")
    end

    @testset "User interface" begin
        include("interface.jl")
    end

    @testset "Unit interface" begin
        include("unit_interface.jl")
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
