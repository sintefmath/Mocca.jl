"""
    water_vapour_pressure(T) → Pa

Saturation pressure of water at temperature `T` [K], from the Antoine equation
with the constants for 1 to 100 °C:

    log₁₀(p_sat / mmHg) = 8.07131 − 1730.63 / (233.426 + T − 273.15)
"""
water_vapour_pressure(T) = 133.322368 * exp10(8.07131 - 1730.63 / (233.426 + T - 273.15))

"""
    TothGAB{N, T} <: AbstractIsotherm

Co-adsorption of CO₂ and water on an amine sorbent (Stampi-Bombelli et al.,
2020, Adsorption 26, 1183): water on a Guggenheim–Anderson–de Boer (GAB)
isotherm in the relative humidity `x = p_w / p_sat(T)`,

    q_w = c_m·c_G·K·x / ((1 − K·x)·(1 + (c_G − 1)·K·x)),

and CO₂ on a Toth isotherm whose capacity and affinity grow with the water
loading,

    q_c = n_s·b·p / (1 + (b·p)^t)^(1/t)
    n_s = n_s0·exp(χ·(1 − T/T0)) / (1 − γ·q_w)
    b   = b0·exp(−ΔH0/(R·T0)·(T0/T − 1))·(1 + β·q_w)
    t   = t0 + α·(1 − T0/T).

Water is taken up independently of CO₂. Other components do not adsorb. The
water loading in `n_s` and `b` is the equilibrium one, with the gas at the
local humidity. `p_sat` is [`water_vapour_pressure`](@ref). The exponent `t`
is kept above 0.05: with `α < 0` it falls with temperature, and at zero the
isotherm breaks down.

Loadings are in mol per m³ of adsorbent, `b0` in 1/Pa and `γ`, `β` in m³/mol.
The heat of adsorption `ΔH` of each component [J/mol, negative] is constant.
Build one by component name with `toth_gab` in the user interface.
"""
# Smallest Toth exponent used; see compute_equilibrium
const TOTH_MIN_EXPONENT = 0.05

struct TothGAB{N, T} <: AbstractIsotherm
    toth::Int             # component on the Toth isotherm (CO₂)
    water::Int            # water, or 0 if the gas has none
    ns0::T                # Toth capacity at T0 [mol/m³ adsorbent]
    b0::T                 # Toth affinity at T0 [1/Pa]
    t0::T                 # Toth exponent at T0 [-]
    α::T                  # temperature dependence of the exponent [-]
    χ::T                  # temperature dependence of the capacity [-]
    ΔH0::T                # isosteric heat of the affinity [J/mol]
    T0::T                 # reference temperature [K]
    γ::T                  # water enhancement of the capacity [m³/mol]
    β::T                  # water enhancement of the affinity [m³/mol]
    cm::T                 # GAB monolayer capacity [mol/m³ adsorbent]
    cG::T                 # GAB constant c_G [-]
    K::T                  # GAB constant K_ads [-]
    ΔH::SVector{N, T}     # heat of adsorption per component [J/mol]
end

function TothGAB(; toth, water, ns0, b0, t0, α, χ = 0.0, ΔH0, T0, γ, β, cm, cG, K, ΔH)
    N = length(ΔH)
    1 <= toth <= N || throw(ArgumentError("toth = $toth is not a component index (1 to $N)"))
    0 <= water <= N && water != toth || throw(ArgumentError("water = $water must be 0 or a component index other than toth"))
    T = promote_type(typeof.((ns0, b0, t0, α, χ, ΔH0, T0, γ, β, cm, cG, K))..., eltype(ΔH), Float64)
    return TothGAB{N, T}(toth, water, ns0, b0, t0, α, χ, ΔH0, T0, γ, β, cm, cG, K, SVector{N, T}(ΔH))
end

"GAB water loading [mol/m³ adsorbent] at relative humidity `x`."
function gab_loading(iso::TothGAB, x)
    Kx = iso.K * x
    return iso.cm * iso.cG * Kx / ((1 - Kx) * (1 + (iso.cG - 1) * Kx))
end

function compute_equilibrium(iso::TothGAB{N}, concentration, temperature) where N
    R = GAS_CONSTANT
    T = temperature
    if iso.water == 0
        q_w = zero(T * concentration[iso.toth])
    else
        x = concentration[iso.water] * R * T / water_vapour_pressure(T)
        q_w = gab_loading(iso, x)
    end
    p = concentration[iso.toth] * R * T
    ns = iso.ns0 * exp(iso.χ * (1 - T / iso.T0)) / (1 - iso.γ * q_w)
    b = iso.b0 * exp(-iso.ΔH0 / (R * iso.T0) * (iso.T0 / T - 1)) * (1 + iso.β * q_w)
    # With α < 0 the exponent falls with temperature and reaches zero (at
    # 399 K for the sorbent of the paper), where the isotherm breaks down;
    # Newton iterations may overshoot that far
    t = max(iso.t0 + iso.α * (1 - iso.T0 / T), TOTH_MIN_EXPONENT)
    # Newton iterations may visit negative partial pressures, where (b·p)^t
    # is not real
    bp = max(b * p, 1e-20)
    q_c, q_w = promote(ns * bp / (1 + bp^t)^(1 / t), q_w)
    return SVector{N}(ntuple(i -> i == iso.toth ? q_c : (i == iso.water ? q_w : zero(q_c)), N))
end

compute_enthalpy(iso::TothGAB, concentration, temperature) = iso.ΔH

# At 20 °C the sorbent of the paper holds 0.003 mol/kg of CO₂ in equilibrium
# with 10⁻¹¹ of it in air at 1 bar
minimum_mole_fraction(iso::TothGAB) = 1e-20

function saturation_capacity(iso::TothGAB{N}) where N
    return SVector{N}(ntuple(i -> i == iso.toth ? iso.ns0 : (i == iso.water ? iso.cm : zero(iso.ns0)), N))
end
