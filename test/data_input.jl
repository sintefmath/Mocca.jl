@testset "Parse JSON" begin
    # Test parsing data from JSON files and dictionaries
    @testset "Parse simple JSON input" begin
        filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input_simple.json")
        (constants_simple_JSON, info_simple_JSON) = Mocca.parse_input(filepath)

        @test constants_simple_JSON isa Mocca.adsorptionConstants{Float64}
        @test info_simple_JSON isa Mocca.processInfo
        @test constants_simple_JSON.molecularMassOfCO2 == 0.04401
        @test constants_simple_JSON.y_feed == [0.15, 0.85]
        @test info_simple_JSON.stage_types == ["adsorption"]
    end

    @testset "Parse detailed JSON input" begin
        filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input.json")
        (constants_JSON, info_JSON) = Mocca.parse_input(filepath)

        @test constants_JSON isa Mocca.adsorptionConstants{Float64}
        @test info_JSON isa Mocca.processInfo
        @test constants_JSON.molecularMassOfCO2 == 0.04401
        @test constants_JSON.y_feed == [0.15, 0.85]
        @test info_JSON.stage_types == ["adsorption"]
    end
    @testset "Parse dictionary input" begin
        (constants_dict, info_dict) = Mocca.parse_input(haghpanah_DCB_input())
        @test constants_dict isa Mocca.adsorptionConstants{Float64}

        @test info_dict isa Mocca.processInfo
        @test constants_dict.molecularMassOfCO2 == 0.04401
        @test constants_dict.y_feed == [0.15, 0.85]
        @test info_dict.stage_types == ["adsorption"]
    end
end