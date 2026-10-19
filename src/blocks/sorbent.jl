# Sorbent block: adsorbed loading, its rate equation, and the terms it adds to
# the gas mass and energy balances. Shared by every SorbentBed.
#
# A SorbentBed system must have `isotherm`, `mass_transfer` and
# `heat_capacity_adsorbed` fields.

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

struct SpecificHeatCapacityAdsorbent <: Jutul.ScalarVariable end

Jutul.@jutul_secondary function update_heat_capacity_adsorbent!(cpa, tv::SpecificHeatCapacityAdsorbent, model::SorbentBedModel, y, ix)
    sys = model.system
    C_pa = sys.heat_capacity_adsorbed
    T = eltype(cpa)
    for cell in ix
        cpa_i = zero(T)
        for component in 1:number_of_components(sys)
            cpa_i += y[component, cell] * C_pa[component]
        end
        cpa[cell] = cpa_i
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

@inline function mass_source(::SorptionUptake, model, state, state0, cell, i, Δt)
    ∂q∂t = Jutul.accumulation_term(state.AdsorbedLoading, state0.AdsorbedLoading, Δt, i, cell)
    return -state.SolidVolume[cell] * ∂q∂t
end

# Components included in the sorption energy terms. Versions up to 0.1.0 only
# included the first component (a loop bound taken from the length-1 energy
# equation buffer); `all_components = false` reproduces that.
@inline _sorption_energy_components(AC, all_components::Bool) = all_components ? axes(AC, 1) : (1:1)

"""
    HeatOfSorption(; all_components = false)

Heat released by sorption, plus the sensible heat carried by the adsorbed
phase: removes `V_s·Σ_i (C_pa·M̄·T + ΔH_i)·∂q_i/∂t` from the gas energy balance.

With `all_components = false` the sum only includes the first component, which
reproduces results from Mocca 0.1.0 and earlier.
"""
struct HeatOfSorption{A} <: AbstractSourceTerm end
HeatOfSorption(; all_components::Bool = false) = HeatOfSorption{all_components}()

@inline function energy_source(::HeatOfSorption{A}, model, state, state0, cell, Δt) where A
    AC = state.AdsorbedLoading
    AC₀ = state0.AdsorbedLoading
    ΔH = state.ΔH
    C_pa = state.C_pa[cell]
    avm = state.AverageMolarMass[cell]
    T = state.Temperature[cell]
    adsorption_term = zero(T)
    for i in _sorption_energy_components(AC, A)
        ∂q∂t = Jutul.accumulation_term(AC, AC₀, Δt, i, cell)
        adsorption_term += (C_pa * avm * T + ΔH[i, cell]) * ∂q∂t
    end
    return -state.SolidVolume[cell] * adsorption_term
end

"Heat capacity of the sorbent solid: `V_s·ρ_s·C_ps`."
struct SolidHeatCapacity <: AbstractSourceTerm end

@inline function heat_capacity(::SolidHeatCapacity, model, state, cell)
    return state.SolidVolume[cell] * state.AdsorbentDensity[1] * state.AdsorbentHeatCapacity[1]
end

"""
    AdsorbedPhaseHeatCapacity(; all_components = false)

Heat capacity of the adsorbed phase: `V_s·C_pa·M̄·Σ_i q_i`. See
[`HeatOfSorption`](@ref) for `all_components`.
"""
struct AdsorbedPhaseHeatCapacity{A} <: AbstractSourceTerm end
AdsorbedPhaseHeatCapacity(; all_components::Bool = false) = AdsorbedPhaseHeatCapacity{all_components}()

@inline function heat_capacity(::AdsorbedPhaseHeatCapacity{A}, model, state, cell) where A
    AC = state.AdsorbedConcentration
    sq = zero(eltype(AC))
    for i in _sorption_energy_components(AC, A)
        sq += AC[i, cell]
    end
    return state.SolidVolume[cell] * state.C_pa[cell] * state.AverageMolarMass[cell] * sq
end
