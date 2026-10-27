"""
    Virial{N, T} <: AbstractIsotherm

Multicomponent Virial isotherm truncated after the second coefficient (Taqvi
and LeVan 1997; Grande et al. 2008), as used for activated carbon by Mondino et
al. (2022, eqs. 12–14 and A.5–A.6). The loading `q*` [mol/kg] is given
implicitly by

    Pᵢ = q*ᵢ/K_H,ᵢ · exp(Σⱼ Aᵢⱼ q*ⱼ + Σⱼₖ Bᵢⱼₖ q*ⱼ q*ₖ)

with `K_H,ᵢ = K∞ᵢ·exp(ΔHᵢ/(R·T))`, `Aᵢ = A0ᵢ + A1ᵢ/T`, `Bᵢ = B0ᵢ + B1ᵢ/T`, and
the mixing rules `Aᵢⱼ = (Aᵢ + Aⱼ)/2`, `Bᵢⱼₖ = (Bᵢ + Bⱼ + Bₖ)/3`.

`compute_equilibrium` solves this by Newton's method on `ln q*` and returns the
loading per particle volume [mol/m³], `density·q*`, as the rest of Mocca uses.
The heat of adsorption is the zero-coverage value `−ΔH`, as in the paper.
"""
struct Virial{N, T} <: AbstractIsotherm
    K∞::SVector{N, T}     # Henry constant at infinite temperature [mol/(kg·Pa)]
    ΔH::SVector{N, T}     # Heat of adsorption at zero coverage, positive [J/mol]
    A0::SVector{N, T}     # First Virial coefficient [kg/mol]
    A1::SVector{N, T}     # Its temperature dependence [K·kg/mol]
    B0::SVector{N, T}     # Second Virial coefficient [kg²/mol²]
    B1::SVector{N, T}     # Its temperature dependence [K·kg²/mol²]
    density::T            # Particle density converting mol/kg to mol/m³ [kg/m³]
end

"""
    Virial(; K∞, ΔH, A0, A1, B0, B1, density)

Construct a [`Virial`](@ref) isotherm in SI units. Mondino et al. (2022) give
`K∞` in mol/(kg·kPa), which is `1e-3` times the value in mol/(kg·Pa).
"""
function Virial(; K∞, ΔH, A0, A1, B0, B1, density)
    N = length(K∞)
    T = promote_type(eltype(K∞), eltype(ΔH), eltype(A0), eltype(A1), eltype(B0), eltype(B1), typeof(density))
    v(x) = SVector{N, T}(x)
    return Virial{N, T}(v(K∞), v(ΔH), v(A0), v(A1), v(B0), v(B1), T(density))
end

function compute_equilibrium(iso::Virial{N}, concentration, temperature) where N
    R = GAS_CONSTANT
    T = temperature
    A = iso.A0 .+ iso.A1 ./ T
    B = iso.B0 .+ iso.B1 ./ T
    # Henry's law loading, the starting point, and its log
    qH = SVector{N}(ntuple(i -> iso.K∞[i] * exp(iso.ΔH[i] / (R * T)) * concentration[i] * R * T, N))
    lnqH = log.(qH)
    u = lnqH
    for _ in 1:30
        q = exp.(u)
        S = sum(q)
        Aq = sum(A .* q)
        Bq = sum(B .* q)
        # Residual of ln q − ln q_Henry + g(q) = 0, with g the exponent above
        g = (A .* S .+ Aq) ./ 2 .+ (B .* S^2 .+ 2 * S * Bq) ./ 3
        r = u .- lnqH .+ g
        # ∂g_i/∂q_j = (A_i + A_j)/2 + 2(B_i·S + Bq + S·B_j)/3, times q_j for ∂/∂u_j
        J = SMatrix{N, N}(ntuple(N * N) do k
            i, j = mod1(k, N), fld1(k, N)
            (i == j) + ((A[i] + A[j]) / 2 + 2 * (B[i] * S + Bq + S * B[j]) / 3) * q[j]
        end)
        du = J \ r
        u = u .- du
        maximum(abs, Jutul.value.(du)) < 1e-13 && break
    end
    return iso.density .* exp.(u)
end

compute_enthalpy(iso::Virial, concentration, temperature) = -iso.ΔH
