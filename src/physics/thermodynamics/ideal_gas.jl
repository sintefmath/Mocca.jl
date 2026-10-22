"""
    IdealGasConstantCp{N, T} <: AbstractGasThermodynamics

Ideal gas whose components have constant molar heat capacities `cp`
[J/(mol·K)].

    IdealGasConstantCp(cp)
    IdealGasConstantCp(; molar_masses, heat_capacities)

The second form takes heat capacities per mass [J/(kg·K)] and the molar masses
[kg/mol] to convert them.
"""
struct IdealGasConstantCp{N, T} <: AbstractGasThermodynamics
    cp::SVector{N, T}    # Molar heat capacity per component [J/(mol·K)]
end

IdealGasConstantCp(cp::AbstractVector) = IdealGasConstantCp(SVector{length(cp)}(float.(cp)))
IdealGasConstantCp(; molar_masses, heat_capacities) = IdealGasConstantCp(collect(molar_masses) .* collect(heat_capacities))

component_heat_capacities(model::IdealGasConstantCp, T) = model.cp
component_enthalpies(model::IdealGasConstantCp, T) = model.cp .* (T - THERMODYNAMIC_REFERENCE_TEMPERATURE)

"""
    IdealGasShomate{N, T} <: AbstractGasThermodynamics

Ideal gas whose components have temperature-dependent heat capacities given by
the Shomate equation, as tabulated in the NIST Chemistry WebBook:

    cp = A + B·t + C·t² + D·t³ + E/t²   [J/(mol·K)],   t = T/(1000 K)

    IdealGasShomate(coefficients)
    IdealGasShomate(; A, B, C, D, E)

`coefficients` holds one `(A, B, C, D, E)` per component; the keyword form one
vector per coefficient. Each set of coefficients is fitted over a range of
temperature, which should cover that of the process.
"""
struct IdealGasShomate{N, T} <: AbstractGasThermodynamics
    A::SVector{N, T}
    B::SVector{N, T}
    C::SVector{N, T}
    D::SVector{N, T}
    E::SVector{N, T}
end

function IdealGasShomate(; A, B, C, D, E)
    N = length(A)
    all(length(x) == N for x in (B, C, D, E)) || throw(ArgumentError("Shomate coefficients A to E must have one value per component"))
    v(x) = SVector{N}(float.(x))
    return IdealGasShomate(v(A), v(B), v(C), v(D), v(E))
end

function IdealGasShomate(coefficients::AbstractVector)
    all(length(c) == 5 for c in coefficients) || throw(ArgumentError("each component needs the five Shomate coefficients (A, B, C, D, E)"))
    coef(k) = [c[k] for c in coefficients]
    return IdealGasShomate(; A = coef(1), B = coef(2), C = coef(3), D = coef(4), E = coef(5))
end

"The Shomate coefficients `(A, B, C, D, E)` of each component of `model`."
shomate_coefficients(model::IdealGasShomate) =
    [(model.A[i], model.B[i], model.C[i], model.D[i], model.E[i]) for i in eachindex(model.A)]

function component_heat_capacities(model::IdealGasShomate, T)
    t = T / 1000
    return model.A .+ model.B .* t .+ model.C .* t^2 .+ model.D .* t^3 .+ model.E ./ t^2
end

# ∫ cp dT from 0 to T, without the integration constant: the enthalpy of the
# Shomate equation without its F and H
function _shomate_enthalpy(model::IdealGasShomate, T)
    t = T / 1000
    return 1000 .* (model.A .* t .+ model.B .* (t^2 / 2) .+ model.C .* (t^3 / 3) .+ model.D .* (t^4 / 4) .- model.E ./ t)
end

component_enthalpies(model::IdealGasShomate, T) =
    _shomate_enthalpy(model, T) .- _shomate_enthalpy(model, THERMODYNAMIC_REFERENCE_TEMPERATURE)
