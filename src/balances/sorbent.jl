# Sorbent block: adsorbed loading, its rate equation, and the terms it adds to
# the gas mass and energy balances. Shared by every SorbentBed.
#
# A SorbentBed system must have `isotherm`, `mass_transfer`, `thermodynamics`
# and `adsorbed_thermodynamics` fields; the last is an AbstractGasThermodynamics
# whose component heat capacities are those of the adsorbed phase.

struct AdsorbedConcentration <: Jutul.VectorVariables
end

function Jutul.minimum_value(::AdsorbedConcentration)
    return 1e-10
end

function Jutul.degrees_of_freedom_per_entity(model::SorbentBedModel, ::AdsorbedConcentration)
    number_of_components(model.system)
end
Jutul.values_per_entity(model::SorbentBedModel, ::AdsorbedConcentration) = number_of_components(model.system)

# The loading as a secondary variable. Every time difference of the loading uses
# this copy rather than the primary variable: when Jutul computes sensitivities
# with respect to parameters it turns the primary variables into parameters,
# which it stores once for both time levels, so a difference of the primary
# variable itself would vanish. Secondary variables are evaluated separately for
# each time level. ColumnConservedEnergy does the same for temperature.
struct AdsorbedLoading <: ComponentVariable end

Jutul.@jutul_secondary function update_adsorbed_loading!(loading, tv::AdsorbedLoading, model::SorbentBedModel, AdsorbedConcentration, ix)
    for cell in ix
        for i in axes(loading, 1)
            loading[i, cell] = AdsorbedConcentration[i, cell]
        end
    end
end

# Rate of change of loading from the mass transfer model [mol/(m³·s)]
struct AdsorptionMassTransfer <: ComponentVariable end

Jutul.@jutul_secondary function update_adsorption_mass_transfer!(
    adsorption_mass_transfer,
    tv::AdsorptionMassTransfer,
    model::SorbentBedModel,
    MolarConcentration,
    Temperature,
    AdsorbedConcentration,
    ix
)
    sys = model.system
    iso = sys.isotherm
    mt = sys.mass_transfer
    N = number_of_components(sys)
    T = eltype(adsorption_mass_transfer)
    for cell in ix
        C = SVector{N, T}(@view MolarConcentration[:, cell])
        q = SVector{N, T}(@view AdsorbedConcentration[:, cell])
        qstar = compute_equilibrium(iso, C, Temperature[cell])
        rate = compute_mass_transfer_rate(mt, C, q, qstar)
        for i in 1:N
            adsorption_mass_transfer[i, cell] = rate[i]
        end
    end
end

"""
The rate of change of the loading as `AdsorptionMassTransfer`, with the
molecular diffusivity of each component following temperature and composition
by the system's `transport`, at the pressure of its reference (the feed). A
bed whose transport properties change uses it in its place.
"""
struct ScaledAdsorptionMassTransfer <: ComponentVariable end

Jutul.@jutul_secondary function update_scaled_adsorption_mass_transfer!(
    adsorption_mass_transfer,
    tv::ScaledAdsorptionMassTransfer,
    model::SorbentBedModel,
    MolarConcentration,
    Temperature,
    AdsorbedConcentration,
    ix
)
    sys = model.system
    iso = sys.isotherm
    mt = sys.mass_transfer
    N = number_of_components(sys)
    T = eltype(adsorption_mass_transfer)
    for cell in ix
        C = SVector{N, T}(@view MolarConcentration[:, cell])
        q = SVector{N, T}(@view AdsorbedConcentration[:, cell])
        qstar = compute_equilibrium(iso, C, Temperature[cell])
        y = C ./ sum(C)
        s = diffusivity_factors(sys.transport, y, Temperature[cell], sys.transport.pressure)
        rate = compute_mass_transfer_rate(mt, C, q, qstar, s)
        for i in 1:N
            adsorption_mass_transfer[i, cell] = rate[i]
        end
    end
end

# Isosteric heat of adsorption per component [J/mol]
struct EnthalpyChange <: ComponentVariable end

Jutul.@jutul_secondary function update_enthalpy_change!(ΔH, tv::EnthalpyChange, model::SorbentBedModel, MolarConcentration, Temperature, ix)
    iso = model.system.isotherm
    N = number_of_components(model.system)
    T = eltype(ΔH)
    for cell in ix
        C = SVector{N, T}(@view MolarConcentration[:, cell])
        ΔH_values = compute_enthalpy(iso, C, Temperature[cell])
        for i in 1:N
            ΔH[i, cell] = ΔH_values[i]
        end
    end
end

"""
Molar heat capacity [J/(mol·K)] of the adsorbed phase in a cell: the average of
the component heat capacities of the system's `adsorbed_thermodynamics`,
weighted by the loadings.
"""
struct AdsorbedMolarHeatCapacity <: Jutul.ScalarVariable end

Jutul.@jutul_secondary function update_heat_capacity_adsorbent!(cpa, tv::AdsorbedMolarHeatCapacity, model::SorbentBedModel, AdsorbedConcentration, Temperature, ix)
    thermo = model.system.adsorbed_thermodynamics
    for cell in ix
        q = view(AdsorbedConcentration, :, cell)
        cpa[cell] = molar_heat_capacity(thermo, q, Temperature[cell]) / sum(q)
    end
end

# Rate equation for the loading: V_s ∂q/∂t − V_s·rate = 0
function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:AdsorbedLoading},
    model::SorbentBedModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    Δt = _scaled_dt(state, Δt)
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    solid_volume = state[:SolidVolume]

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] =
            solid_volume[self_cell] * ∂M∂t -
            solid_volume[self_cell] * forcing_term[component, self_cell]
    end
end

"Gas taken up by the sorbent: removes `V_s·∂q_i/∂t` from the gas mass balance."
struct SorptionUptake <: AbstractSourceTerm end

@inline function Framework.mass_source(::SorptionUptake, model, state, state0, cell, i, Δt)
    ∂q∂t = Jutul.accumulation_term(state.AdsorbedLoading, state0.AdsorbedLoading, Δt, i, cell)
    return -state.SolidVolume[cell] * ∂q∂t
end

"""
Heat released by sorption, and the sensible heat that the moles taken up carry
from the gas to the adsorbed phase: removes
`V_s·Σ_i (C_pg·T + Δh_i + ΔH_i)·∂q_i/∂t` from the gas energy balance.

The gas energy terms of the column, `(C_pg/R)·(V_f·∂P/∂t + ∇·(vP))`, count
`C_pg·T` for every mole the gas loses, which `C_pg·T` here gives back, so that
only the heat of sorption `ΔH_i` and `Δh_i`, the enthalpy of the adsorbed
component less that of the gas, are left. `Δh_i = (c_pa,i − c_pg,i)·T` for
constant heat capacities (Haghpanah et al. 2013 write `c_pa = c_pg`), and is
zero when the adsorbed phase has the heat capacities of the gas.
"""
struct HeatOfSorption <: AbstractSourceTerm end

@inline function Framework.energy_source(::HeatOfSorption, model, state, state0, cell, Δt)
    AC = state.AdsorbedLoading
    AC₀ = state0.AdsorbedLoading
    ΔH = state.ΔH
    T = state.Temperature[cell]
    gas = state.C_pg[cell] * T
    Δh = _adsorbed_excess_enthalpies(model.system, T)
    adsorption_term = zero(T)
    for i in axes(AC, 1)
        ∂q∂t = Jutul.accumulation_term(AC, AC₀, Δt, i, cell)
        adsorption_term += (gas + Δh[i] + ΔH[i, cell]) * ∂q∂t
    end
    return -state.SolidVolume[cell] * adsorption_term
end

# Enthalpy of each adsorbed component less that of the gas component [J/mol]:
# (c_pa,i − c_pg,i)·T for constant heat capacities, as in Haghpanah et al.;
# otherwise that difference at the reference temperature, carried to T by the
# integrals of the heat capacities (Kirchhoff's law). Exactly zero when the
# two models are the same.
function _adsorbed_excess_enthalpies(sys, T)
    gas, ads = sys.thermodynamics, sys.adsorbed_thermodynamics
    T0 = THERMODYNAMIC_REFERENCE_TEMPERATURE
    return (component_heat_capacities(ads, T0) .- component_heat_capacities(gas, T0)) .* T0 .+
        component_enthalpies(ads, T) .- component_enthalpies(gas, T)
end

"Heat capacity of the sorbent solid: `V_s·ρ_s·C_ps`."
struct SolidHeatCapacity <: AbstractSourceTerm end

@inline function Framework.heat_capacity(::SolidHeatCapacity, model, state, cell)
    return state.SolidVolume[cell] * state.AdsorbentDensity[1] * state.AdsorbentHeatCapacity[1]
end

"Heat capacity of the adsorbed phase: `V_s·C_pa·Σ_i q_i = V_s·Σ_i c_pa,i·q_i`."
struct AdsorbedPhaseHeatCapacity <: AbstractSourceTerm end

@inline function Framework.heat_capacity(::AdsorbedPhaseHeatCapacity, model, state, cell)
    AC = state.AdsorbedConcentration
    sq = zero(eltype(AC))
    for i in axes(AC, 1)
        sq += AC[i, cell]
    end
    return state.SolidVolume[cell] * state.C_pa[cell] * sq
end
