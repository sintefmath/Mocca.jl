using LinearAlgebra
using Tullio
function Jutul.convergence_criterion(model::Jutul.SimulationModel{D,S}, storage, eq::Jutul.ConservationLaw, eq_s, r; dt=1.0, update_report = missing) where {D,S<:AdsorptionSystem}
    n = Jutul.number_of_equations_per_entity(model, eq)
    @tullio max e[i] := abs(r[i, j])
    if n == 1
        names = "R"
    else
        names = map(i -> "R_$i", 1:n)
    end
    R = (AbsMax = (errors = e, names = names), )
    return R
end
