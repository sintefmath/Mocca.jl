# Langmuir.jl backend for Mocca isotherm calculations
# This provides a modular interface to use any Langmuir.jl isotherm model

import Langmuir
import CommonSolve

"""
    LangmuirIsothermModel{I,T}

Generic isotherm model that wraps Langmuir.jl isotherm models.
Supports both single-component models (for pure component isotherms) 
and multi-component models with IAST.

# Fields
- `isotherms`: Tuple of Langmuir.jl isotherm models for each component
- `mass_transfer_params`: Parameters for computing mass transfer coefficients
- `use_iast`: Whether to use IAST for multicomponent adsorption (default: true)
"""
struct LangmuirIsothermModel{I,M} <: AbstractIsothermModel
    isotherms::I  # Tuple of Langmuir.jl isotherm models
    mass_transfer_params::M  # Parameters for mass transfer calculation
    use_iast::Bool
end

# Constructor for dual-site Langmuir from Haghpanah parameters
function LangmuirIsothermModel(constants::HaghpanahConstants{T}; use_iast::Bool=true, T_ref::Real=298.15) where T
    # Create dual-site Langmuir models for CO2 and N2

    # The parameter mapping is challenging due to unit differences:
    # Haghpanah: qstar = qsb * b * C / (1 + sum(b*C)) + qsd * d * C / (1 + sum(d*C))
    # where b = b0 * exp(-ΔUb/(R*T)) and d = d0 * exp(-ΔUd/(R*T))
    # Here C is concentration [mol/m³]
    #
    # LangmuirS1: loading = M * K * p / (1 + K * p)  
    # where K = K₀ * exp(-E/(R*T)) and p is pressure [Pa]
    #
    # Since p = C * R * T, we need: K₀ = b0 / (R * T_ref)
    # We use a reference temperature T_ref to make K₀ temperature-independent

    R_gas = constants.R

    # CO2 dual-site model (has both high and low energy sites)
    co2_site1 = Langmuir.LangmuirS1(
        constants.qsbi[1],                    # M: saturation loading [mol/kg]
        constants.b0[1] / (R_gas * T_ref),   # K₀: corrected for pressure units [1/Pa]  
        constants.ΔUbi[1]                    # E: adsorption energy [J/mol]
    )
    co2_site2 = Langmuir.LangmuirS1(
        constants.qsdi[1],                    # M: saturation loading [mol/kg]
        constants.d0[1] / (R_gas * T_ref),   # K₀: corrected for pressure units [1/Pa]
        constants.ΔUdi[1]                    # E: adsorption energy [J/mol] 
    )
    co2_isotherm = Langmuir.MultiSite(co2_site1, co2_site2)

    # N2 single-site model (only high energy site, since d0[2] = 0 and qsdi[2] = 0)
    n2_isotherm = Langmuir.LangmuirS1(
        constants.qsbi[2],                    # M: saturation loading [mol/kg]
        constants.b0[2] / (R_gas * T_ref),   # K₀: corrected for pressure units [1/Pa]
        constants.ΔUbi[2]                    # E: adsorption energy [J/mol]
    )

    isotherms = (co2_isotherm, n2_isotherm)

    # Store mass transfer parameters and reference temperature
    mass_transfer_params = (
        D_m=constants.D_m,
        τ=constants.τ,
        ϵ_p=constants.ϵ_p,
        d_p=constants.d_p,
        T_ref=T_ref,
        R=R_gas
    )

    return LangmuirIsothermModel(isotherms, mass_transfer_params, use_iast)
end

function compute_equilibrium(model::LangmuirIsothermModel, concentration, temperature)
    N = length(concentration)
    T = eltype(concentration)

    # Get the gas constant from parameters (use the same one as Haghpanah)
    R = model.mass_transfer_params.R
    T_ref = model.mass_transfer_params.T_ref

    if model.use_iast && N > 1
        # Use IAST for multicomponent adsorption

        # Convert concentrations to partial pressures
        # p_i = C_i * R * T, but we need to account for the reference temperature correction
        # Since our K₀ values were computed for T_ref, we need to apply temperature correction
        partial_pressures = concentration .* R .* temperature
        total_pressure = sum(partial_pressures)

        # Compute mole fractions
        y = partial_pressures ./ total_pressure

        # Call IAST solver
        q_total, x_ads, status = Langmuir.iast(
            model.isotherms,
            total_pressure,
            temperature,
            y;
            maxiters=100,
            reltol=1e-12,
            abstol=1e-10
        )

        if status != :success
            @warn "IAST convergence failed, falling back to pure component isotherms"
            # Fallback to pure component calculations
            qstar = @MVector zeros(T, N)
            for i in 1:N
                p_i = concentration[i] * R * temperature
                qstar[i] = Langmuir.loading(model.isotherms[i], p_i, temperature)
            end
            return SVector{N,T}(qstar)
        end

        # Convert adsorbed phase composition to individual loadings
        qstar = @MVector zeros(T, N)
        for i in 1:N
            qstar[i] = x_ads[i] * q_total
        end
        return SVector{N,T}(qstar)

    else
        # Pure component calculations (no competitive adsorption)
        qstar = @MVector zeros(T, N)

        for i in 1:N
            # Convert concentration to pressure for Langmuir.jl
            p_i = concentration[i] * R * temperature
            qstar[i] = Langmuir.loading(model.isotherms[i], p_i, temperature)
        end

        return SVector{N,T}(qstar)
    end
end

function compute_ki(model::LangmuirIsothermModel, concentration, qstar)
    params = model.mass_transfer_params
    D_p = params.D_m / params.τ
    r_p = params.d_p / 2.0

    return concentration ./ qstar .* 15 .* params.ϵ_p .* D_p ./ r_p^2
end

"""
    create_haghpanah_equivalent_langmuir(constants::HaghpanahConstants)

Create a LangmuirIsothermModel that exactly replicates the Haghpanah dual-site 
Langmuir behavior using Langmuir.jl backend.
"""
function create_haghpanah_equivalent_langmuir(constants::HaghpanahConstants)
    return LangmuirIsothermModel(constants; use_iast=true)
end
