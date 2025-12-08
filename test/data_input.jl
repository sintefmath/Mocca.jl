@testset "Parse simple JSON" begin
    # Test parsing data from JSON files and dictionaries
    filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input_simple.json")
    input_pars_simple_JSON = Mocca.parse_input(filepath)

    filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input.json")
    input_pars_JSON = Mocca.parse_input(filepath)

    input_pars_dict = Mocca.parse_input(haghpanah_DCB_input())

    @test input_pars_simple_JSON isa Mocca.adsorptionInput{Float64}
    @test input_pars_JSON isa Mocca.adsorptionInput{Float64}
    @test input_pars_dict isa Mocca.adsorptionInput{Float64}

end