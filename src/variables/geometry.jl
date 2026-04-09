struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    dx = data_domain[:dx]
    if T == :in
        r = first(data_domain[:r_in, Column()])
    else
        @assert T == :out
        r = first(data_domain[:r_out, Column()])
    end
    return 2π .* r .* dx  # Lateral surface area of cylindrical shell
end

struct CellDx <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::CellDx, symb)
    return copy(data_domain[:dx])
end
