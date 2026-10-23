"""
    SutherlandTransport{N, T} <: AbstractGasTransport

Transport properties of a mixture of ideal gases from component data:

- viscosity of each component by Sutherland's law,
  `μᵢ = μ₀ᵢ·(T/T₀)^(3/2)·(T₀ + Sᵢ)/(T + Sᵢ)`, mixed by Herning and Zipperer,
  `μ = Σᵢ yᵢ μᵢ √Mᵢ / Σᵢ yᵢ √Mᵢ`;
- thermal conductivity of each component by a power law,
  `kᵢ = k₀ᵢ·(T/T₀)^nᵢ`, mixed by Wassiljewa's equation with the coefficients of
  Mason and Saxena, `k = Σᵢ yᵢ kᵢ / Σⱼ yⱼ Φᵢⱼ`, with
  `Φᵢⱼ = [1 + (μᵢ/μⱼ)^(1/2) (Mⱼ/Mᵢ)^(1/4)]² / [8 (1 + Mᵢ/Mⱼ)]^(1/2)`;
- diffusivity in the pores from molecular and Knudsen diffusion in series
  (Bosanquet), `1/D ∝ P·(T₀/T)^m + P_K·(T₀/T)^(1/2)`: molecular diffusion goes
  as `T^m/P`, with `m = 1.75` as in the correlation of Fuller, Schettler and
  Giddings, and Knudsen diffusion as `T^(1/2)`, whatever the pressure. `P_K`
  [Pa] is the pressure at which the two are equal at `T₀`, where the mean free
  path of the gas is about the pore diameter (some 0.1 bar for pores of
  0.5 μm). With `P_K = 0`, the default, the diffusion is molecular.

`T₀` is [`THERMODYNAMIC_REFERENCE_TEMPERATURE`](@ref) (298.15 K).

    SutherlandTransport(; viscosity, sutherland, thermal_conductivity,
        conductivity_exponent, molar_masses, diffusivity_exponent = 1.75,
        knudsen_pressure = 0.0)

with one value per component of each: viscosity [Pa·s] and thermal
conductivity [W/(m·K)] at 298.15 K, Sutherland constant [K], conductivity
exponent and molar mass [kg/mol].
"""
struct SutherlandTransport{N, T} <: AbstractGasTransport
    μ₀::SVector{N, T}     # Viscosity at T₀ [Pa·s]
    S::SVector{N, T}      # Sutherland constant [K]
    k₀::SVector{N, T}     # Thermal conductivity at T₀ [W/(m·K)]
    n::SVector{N, T}      # Exponent of the thermal conductivity
    M::SVector{N, T}      # Molar mass [kg/mol]
    m::T                  # Exponent of the temperature in molecular diffusion
    P_K::T                # Pressure at which Knudsen and molecular diffusion are equal at T₀ [Pa]
end

function SutherlandTransport(; viscosity, sutherland, thermal_conductivity, conductivity_exponent, molar_masses,
        diffusivity_exponent = 1.75, knudsen_pressure = 0.0)
    knudsen_pressure >= 0 || throw(ArgumentError("knudsen_pressure must not be negative"))
    N = length(viscosity)
    all(length(x) == N for x in (sutherland, thermal_conductivity, conductivity_exponent, molar_masses)) ||
        throw(ArgumentError("SutherlandTransport needs one value per component of each property"))
    v(x) = SVector{N}(float.(x))
    return SutherlandTransport(v(viscosity), v(sutherland), v(thermal_conductivity), v(conductivity_exponent), v(molar_masses),
        float(diffusivity_exponent), float(knudsen_pressure))
end

"Viscosity [Pa·s] of each component of `model` at temperature `T` [K]."
function component_viscosities(model::SutherlandTransport, T)
    T0 = THERMODYNAMIC_REFERENCE_TEMPERATURE
    return model.μ₀ .* (T / T0)^1.5 .* (T0 .+ model.S) ./ (T .+ model.S)
end

"Thermal conductivity [W/(m·K)] of each component of `model` at temperature `T` [K]."
component_thermal_conductivities(model::SutherlandTransport, T) =
    model.k₀ .* (T / THERMODYNAMIC_REFERENCE_TEMPERATURE) .^ model.n

function mixture_viscosity(model::SutherlandTransport, y, T)
    μ = component_viscosities(model, T)
    num = zero(y[1] * μ[1])
    den = zero(y[1] * one(eltype(model.M)))
    for i in eachindex(μ)
        w = y[i] * sqrt(model.M[i])
        num += w * μ[i]
        den += w
    end
    return num / den
end

function mixture_thermal_conductivity(model::SutherlandTransport, y, T)
    μ = component_viscosities(model, T)
    k = component_thermal_conductivities(model, T)
    M = model.M
    s = zero(y[1] * k[1])
    for i in eachindex(k)
        d = zero(y[1] * μ[1] / μ[1])
        for j in eachindex(k)
            Φ = (1 + sqrt(μ[i] / μ[j]) * sqrt(sqrt(M[j] / M[i])))^2 / sqrt(8 * (1 + M[i] / M[j]))
            d += y[j] * Φ
        end
        s += y[i] * k[i] / d
    end
    return s
end

function diffusivity_scale(model::SutherlandTransport, T, P)
    θ = THERMODYNAMIC_REFERENCE_TEMPERATURE / T
    return 1 / (P * θ^model.m + model.P_K * sqrt(θ))
end
