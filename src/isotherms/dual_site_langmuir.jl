using StaticArrays

"""
    DualSiteLangmuir{N, T} <: AbstractIsotherm

Dual-site Langmuir isotherm (Haghpanah et al. 2013), concentration basis.

    q*ᵢ = qsbᵢ·bᵢ·Cᵢ/(1 + Σⱼ bⱼ·Cⱼ) + qsdᵢ·dᵢ·Cᵢ/(1 + Σⱼ dⱼ·Cⱼ)

Two adsorption sites per component:
- Site b (high energy): saturation `qsb`, affinity `b(T) = b0·exp(-ΔUb/(R·T))`
- Site d (low energy): saturation `qsd`, affinity `d(T) = d0·exp(-ΔUd/(R·T))`

Enthalpy is state-independent, precomputed as a weighted average of site energies
normalized by the primary adsorbate's total saturation capacity.
"""
struct DualSiteLangmuir{N, T} <: AbstractIsotherm
    # Site b (high energy)
    qsb::SVector{N, T}    # Saturation loading [mol/kg]
    b0::SVector{N, T}     # Pre-exponential factor [m³/mol]
    ΔUb::SVector{N, T}    # Adsorption energy [J/mol]
    # Site d (low energy)
    qsd::SVector{N, T}    # Saturation loading [mol/kg]
    d0::SVector{N, T}     # Pre-exponential factor [m³/mol]
    ΔUd::SVector{N, T}    # Adsorption energy [J/mol]
    # Constants
    R::T                   # Gas constant [J/(mol·K)]
    T0::T                  # Reference temperature [K]
    ΔH::SVector{N, T}      # Isosteric heat of adsorption [J/mol]
end

function DualSiteLangmuir(constants::ConstantsStruct)
    T = typeof(constants.R)
    N = 2
    # Weighted average of site energies, normalized by primary adsorbate saturation capacity
    sumq = constants.qsbi[1] + constants.qsdi[1]
    ΔH = SVector{N}(ntuple(i ->
        (constants.qsbi[i] * (constants.ΔUbi[i] - constants.R * constants.T0) +
         constants.qsdi[i] * (constants.ΔUdi[i] - constants.R * constants.T0)) / sumq,
    N))
    return DualSiteLangmuir{N, T}(
        constants.qsbi,
        constants.b0,
        constants.ΔUbi,
        constants.qsdi,
        constants.d0,
        constants.ΔUdi,
        constants.R,
        constants.T0,
        ΔH
    )
end

function compute_equilibrium(iso::DualSiteLangmuir{N}, concentration, temperature) where N
    # Temperature-dependent affinity constants
    b = SVector{N}(ntuple(i -> iso.b0[i] * exp(-iso.ΔUb[i] / (iso.R * temperature)), N))
    d = SVector{N}(ntuple(i -> iso.d0[i] * exp(-iso.ΔUd[i] / (iso.R * temperature)), N))

    # Competitive adsorption denominators
    bC_sum = zero(b[1] * concentration[1])
    dC_sum = zero(d[1] * concentration[1])
    @inbounds for i in 1:N
        bC_sum += b[i] * concentration[i]
        dC_sum += d[i] * concentration[i]
    end

    # Equilibrium loading per component
    return SVector{N}(ntuple(i ->
        iso.qsb[i] * b[i] * concentration[i] / (one(bC_sum) + bC_sum) +
        iso.qsd[i] * d[i] * concentration[i] / (one(dC_sum) + dC_sum),
    N))
end

compute_enthalpy(iso::DualSiteLangmuir, concentration, temperature) = iso.ΔH
