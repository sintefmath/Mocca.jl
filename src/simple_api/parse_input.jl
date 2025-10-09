using Mocca

function parse_input(input_dict::Dict{String, Any})
        return Mocca.setup_case_from_dict(input_dict)
end

function parse_input(filepath::String)
    using JSON3
    input_dict = JSON3.read(open(filepath), Dict{String, Any})
    return Mocca.setup_case_from_dict(input_dict)
end