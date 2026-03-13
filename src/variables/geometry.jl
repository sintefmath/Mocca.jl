struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    T::Symbol
    dx = data_domain[:dx]
    sys = model.system
    if T == :in
        F = x -> area_wall_in(sys, x)
    else
        @assert T == :out
        F = x -> area_wall_out(sys, x)
    end
    return F.(dx)
end

struct CellDx <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::CellDx, symb)
    return copy(data_domain[:dx])
end
