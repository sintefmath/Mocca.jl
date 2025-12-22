@testset "Parse simple JSON" begin
    # Test parsing data from JSON files and dictionaries
    filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input_simple.json")
    (constants_simple_JSON, info_simple_JSON) = Mocca.parse_input(filepath)

    filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input.json")
    (constants_JSON, info_JSON) = Mocca.parse_input(filepath)

    (constants_dict, info_dict) = Mocca.parse_input(haghpanah_DCB_input())

    @test constants_simple_JSON isa Mocca.adsorptionConstants{Float64}
    @test constants_JSON isa Mocca.adsorptionConstants{Float64}
    @test constants_dict isa Mocca.adsorptionConstants{Float64}

    @test info_simple_JSON isa Mocca.processInfo
    @test info_JSON isa Mocca.processInfo
    @test info_dict isa Mocca.processInfo
end