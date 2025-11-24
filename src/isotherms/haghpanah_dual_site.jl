# Original Haghpanah dual-site Langmuir isotherm implementation
# This preserves the exact behavior from the original Mocca implementation

"""
    HaghpanahDualSiteIsotherm{T}

Dual-site Langmuir isotherm model using Haghpanah et al. 2013 parameterization.
This exactly replicates the original compute_equilibrium function behavior.
"""
struct HaghpanahDualSiteIsotherm{T<:ConstantsStruct} <: AbstractIsothermModel
    p::T  # Parameters struct (e.g., HaghpanahConstants)
end

function compute_equilibrium(model::HaghpanahDualSiteIsotherm, concentration, temperature)
    sys_p = model.p
    N = length(concentration)
    T = eltype(concentration)
    qstar = @MVector zeros(T, N)
    b = @MVector zeros(T, N)
    d = @MVector zeros(T, N)
    bC_sum = zero(T)
    dC_sum = zero(T)

    @inbounds for i in 1:N
        b_i = sys_p.b0[i] * exp(-sys_p.ΔUbi[i] / (sys_p.R * temperature))
        d_i = sys_p.d0[i] * exp(-sys_p.ΔUdi[i] / (sys_p.R * temperature))

        b[i] = b_i
        d[i] = d_i

        bC_sum += b_i * concentration[i]
        dC_sum += d_i * concentration[i]
    end

    @inbounds for i in 1:N
        qstar[i] = (sys_p.qsbi[i] * b[i] * concentration[i] / (one(T) + bC_sum) +
                   sys_p.qsdi[i] * d[i] * concentration[i] / (one(T) + dC_sum))
    end
    return SVector{N, T}(qstar)
end

function compute_ki(model::HaghpanahDualSiteIsotherm, concentration, qstar)
    sys_p = model.p
    D_p = sys_p.D_m / sys_p.τ
    r_p = sys_p.d_p / 2.0

    return concentration ./ qstar .* 15 .* sys_p.ϵ_p .* D_p ./ r_p^2
end
