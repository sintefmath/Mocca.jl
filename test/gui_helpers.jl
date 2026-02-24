@testset "GUI Helper Functions" begin

    @testset "Load default parameters" begin
        params = Mocca.load_default_parameters()
        @test params isa Dict{String, Any}
        @test haskey(params, "physicalConstants")
        @test haskey(params, "columnProps")
        @test haskey(params, "processSpecification")
        @test haskey(params, "simulation")
        @test haskey(params, "solver")
        @test params["columnProps"]["L"] == 1.0
        @test params["processSpecification"]["num_cycles"] == 500
        @test params["processSpecification"]["stage_types"] == ["pressurisation", "adsorption", "blowdown", "evacuation"]
    end

    @testset "Load parameters from file" begin
        # Simple format
        simple_path = joinpath(@__DIR__, "../models/json/haghpanah_cyclic_input_simple.json")
        params_simple = Mocca.load_parameters_from_file(simple_path)
        @test params_simple isa Dict{String, Any}
        @test params_simple["columnProps"]["L"] == 1.0

        # Detailed format
        detailed_path = joinpath(@__DIR__, "../models/json/haghpanah_cyclic_input.json")
        params_detailed = Mocca.load_parameters_from_file(detailed_path)
        @test params_detailed isa Dict{String, Any}
        @test haskey(params_detailed["columnProps"]["L"], "value")
    end

    @testset "Export parameters to file" begin
        params = Mocca.load_default_parameters()
        tmpfile = tempname() * ".json"
        Mocca.export_parameters_to_file(tmpfile, params)
        @test isfile(tmpfile)

        # Reimport and verify
        reimported = Mocca.load_parameters_from_file(tmpfile)
        @test reimported["columnProps"]["L"] == params["columnProps"]["L"]
        @test reimported["feedProps"]["y_feed"] == params["feedProps"]["y_feed"]
        @test reimported["processSpecification"]["stage_types"] == params["processSpecification"]["stage_types"]
        @test reimported["processSpecification"]["num_cycles"] == params["processSpecification"]["num_cycles"]

        rm(tmpfile)
    end

    @testset "value_to_string" begin
        @test Mocca.value_to_string(1.0) == "1.0"
        @test Mocca.value_to_string(500) == "500"
        @test Mocca.value_to_string("test") == "test"
        @test Mocca.value_to_string([0.15, 0.85]) == "0.15, 0.85"
        @test Mocca.value_to_string([8.65e-7, 2.5e-6]) == "8.65e-7, 2.5e-6"
        @test Mocca.value_to_string(Dict()) == "{}"
    end

    @testset "string_to_value" begin
        # Float
        @test Mocca.string_to_value("1.0", 0.0) == 1.0
        @test Mocca.string_to_value("8.65e-7", 0.0) == 8.65e-7

        # Integer
        @test Mocca.string_to_value("500", 0) == 500

        # String
        @test Mocca.string_to_value("test", "") == "test"

        # Float array
        @test Mocca.string_to_value("0.15, 0.85", [0.0, 0.0]) == [0.15, 0.85]
        @test Mocca.string_to_value("8.65e-7, 2.5e-6", [0.0, 0.0]) == [8.65e-7, 2.5e-6]
        @test Mocca.string_to_value("-36641.21, -15800.0", [0.0, 0.0]) == [-36641.21, -15800.0]

        # String array
        @test Mocca.string_to_value("pressurisation, adsorption", ["a", "b"]) == ["pressurisation", "adsorption"]

        # Dict (passthrough)
        @test Mocca.string_to_value("{}", Dict()) == Dict()
    end

    @testset "Export then parse roundtrip" begin
        # Load defaults, export, reimport, parse into simulation structs
        params = Mocca.load_default_parameters()
        tmpfile = tempname() * ".json"
        Mocca.export_parameters_to_file(tmpfile, params)

        reimported = Mocca.load_parameters_from_file(tmpfile)
        constants, info = Mocca.parse_input(reimported)

        @test constants isa Mocca.adsorptionConstants{Float64}
        @test info isa Mocca.processInfo
        @test constants.L == 1.0
        @test constants.Φ == 0.37
        @test constants.y_feed == [0.15, 0.85]
        @test info.num_cycles == 500
        @test info.stage_types == ["pressurisation", "adsorption", "blowdown", "evacuation"]
        @test info.stage_durations == [15.0, 15.0, 30.0, 40.0]
        @test info.system_type == "TwoComponentAdsorptionSystem"

        rm(tmpfile)
    end
end
