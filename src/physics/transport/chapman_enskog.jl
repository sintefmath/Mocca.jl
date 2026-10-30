"""
    ChapmanEnskogTransport{N, T} <: AbstractGasTransport

Transport properties of a mixture of ideal gases by the kinetic theory of
Chapman and Enskog, with the Lennard-Jones parameters of each component and
the collision integrals of Neufeld et al. (Poling, Prausnitz and O'Connell,
*The Properties of Gases and Liquids*, 5th ed., §9-3, §9-5 and §11-3):

- viscosity of each component, `μᵢ = 2.6693e-6·√(Mᵢ T) / (σᵢ² Ω_μ(T/εᵢ))`
  [Pa·s], with `Mᵢ` in g/mol and `σᵢ` in Å, mixed by Wilke,
  `μ = Σᵢ yᵢ μᵢ / Σⱼ yⱼ Φᵢⱼ`, `Φᵢⱼ = [1 + (μᵢ/μⱼ)^(1/2) (Mⱼ/Mᵢ)^(1/4)]² / [8 (1 + Mᵢ/Mⱼ)]^(1/2)`;
- binary diffusivities, `Dᵢⱼ = 1.858e-7·T^(3/2)·√(1/Mᵢ + 1/Mⱼ) / (P σᵢⱼ² Ω_D(T/εᵢⱼ))`
  [m²/s], with `P` in atm, `σᵢⱼ = (σᵢ + σⱼ)/2` and `εᵢⱼ = √(εᵢ εⱼ)`, and the
  diffusivity of each component in the mixture by Wilke,
  `Dᵢ = Σ_{j≠i} yⱼ / Σ_{j≠i} (yⱼ/Dᵢⱼ)`.

    ChapmanEnskogTransport(; sigma, epsilon, molar_masses)

with one value per component of each: collision diameter `sigma` [Å], well
depth `epsilon` as ε/k [K], and molar mass [kg/mol].
"""
struct ChapmanEnskogTransport{N, T} <: AbstractGasTransport
    σ::SVector{N, T}      # Lennard-Jones collision diameter [Å]
    ε::SVector{N, T}      # Lennard-Jones well depth ε/k [K]
    M::SVector{N, T}      # Molar mass [g/mol]
end

function ChapmanEnskogTransport(; sigma, epsilon, molar_masses)
    N = length(sigma)
    length(epsilon) == N && length(molar_masses) == N ||
        throw(ArgumentError("ChapmanEnskogTransport needs one sigma, epsilon and molar mass per component"))
    v(x) = SVector{N}(float.(x))
    return ChapmanEnskogTransport(v(sigma), v(epsilon), 1000 .* v(molar_masses))
end

# Collision integrals of Neufeld, Janzen and Aziz (1972) at the reduced
# temperature T* = kT/ε
@inline _omega_viscosity(Ts) = 1.16145 / Ts^0.14874 + 0.52487 / exp(0.77320 * Ts) + 2.16178 / exp(2.43787 * Ts)
@inline _omega_diffusion(Ts) = 1.06036 / Ts^0.15610 + 0.19300 / exp(0.47635 * Ts) +
    1.03587 / exp(1.52996 * Ts) + 1.76474 / exp(3.89411 * Ts)

"Viscosity [Pa·s] of each component of `model` at temperature `T` [K]."
component_viscosities(model::ChapmanEnskogTransport, T) =
    2.6693e-6 .* sqrt.(model.M .* T) ./ (model.σ .^ 2 .* _omega_viscosity.(T ./ model.ε))

function mixture_viscosity(model::ChapmanEnskogTransport, y, T)
    μ = component_viscosities(model, T)
    M = model.M
    s = zero(y[1] * μ[1])
    for i in eachindex(μ)
        d = zero(y[1] * μ[1] / μ[1])
        for j in eachindex(μ)
            Φ = (1 + sqrt(μ[i] / μ[j]) * sqrt(sqrt(M[j] / M[i])))^2 / sqrt(8 * (1 + M[i] / M[j]))
            d += y[j] * Φ
        end
        s += y[i] * μ[i] / d
    end
    return s
end

"""
    binary_diffusivity(model::ChapmanEnskogTransport, i, j, T, P)

Diffusivity [m²/s] of components `i` and `j` of `model` in each other at
temperature `T` [K] and pressure `P` [Pa].
"""
function binary_diffusivity(model::ChapmanEnskogTransport, i, j, T, P)
    σ = (model.σ[i] + model.σ[j]) / 2
    ε = sqrt(model.ε[i] * model.ε[j])
    return 1.858e-7 * T^1.5 * sqrt(1 / model.M[i] + 1 / model.M[j]) / (P / 101325 * σ^2 * _omega_diffusion(T / ε))
end

# A floor on the mole fractions of the other components. Where a component is
# nearly pure, its mixture diffusivity is a mean of its binary diffusivities
# weighted by the traces of the others, and swings between them as the ratio of
# those traces changes: with three components and traces of 1e-10 its
# derivatives reach 1e10, and Newton's method fails. The floor makes the mean
# of a nearly pure component that of its binary diffusivities, and changes no
# weight by more than 1% for mole fractions above 1e-4. Mole fractions a
# Newton iteration takes below zero count as zero.
const _TRACE = 1e-6

function molecular_diffusivities(model::ChapmanEnskogTransport{N}, y, T, P) where N
    N == 1 && return SVector{1}(binary_diffusivity(model, 1, 1, T, P))
    return SVector{N}(ntuple(Val(N)) do i
        num = zero(y[1])
        den = zero(y[1] / binary_diffusivity(model, 1, 2, T, P))
        for j in 1:N
            j == i && continue
            w = max(y[j], zero(y[j])) + _TRACE
            num += w
            den += w / binary_diffusivity(model, i, j, T, P)
        end
        num / den
    end)
end
