# Sorbent advection: in a unit that `moves_sorbent`, the sorbent moves down
# through the unit, from the top cell (nc) to the bottom cell (1), carrying its
# loading, its heat and the gas in its pores. Upwinding is fixed by the
# direction of motion: what leaves cell c through its bottom face takes the
# value reconstructed from cells c and c + 1 with the system's limiter (see
# limiters.jl), and the sorbent entering the top cell comes from the solid
# inlet. The faces below the top cell and below the bottom cell are first
# order, so that a unit's inlet and outlet carry the values of its end cells,
# as the cross terms between units expect.
#
# Q_s = SorbentMassFlow / AdsorbentDensity is the volume flow of particles
# [m³/s]; the loading is per particle volume, so Q_s·q is a molar flow.
#
# Part of the change in loading leaves with the sorbent, so sorption enters the
# gas balances through the mass transfer rate (SorptionRate, SorptionRateHeat)
# rather than through ∂q/∂t (SorptionUptake, HeatOfSorption). With no sorbent
# flow the two forms agree once the loading equation is solved.

"""
    SorbentMassFlow

Mass flow of sorbent through a unit that moves its sorbent [kg/s], from top to
bottom. Units passing sorbent to each other must have the same value.
"""
struct SorbentMassFlow <: Jutul.ScalarVariable end
Jutul.associated_entity(::SorbentMassFlow) = Unit()
Jutul.minimum_value(::SorbentMassFlow) = 0.0
function Jutul.default_parameter_values(data_domain, model, ::SorbentMassFlow, symb)
    return copy(data_domain[:sorbent_mass_flow, Unit()])
end

"""
    SolidInletLoading

Loading of each component on the sorbent entering the top of a unit that
moves its sorbent [mol/m³ of particles]. Not used when the top is fed by
another unit.
"""
struct SolidInletLoading <: ComponentVariable end
Jutul.associated_entity(::SolidInletLoading) = Unit()
Jutul.minimum_value(::SolidInletLoading) = 0.0
function Jutul.default_parameter_values(data_domain, model, ::SolidInletLoading, symb)
    return copy(data_domain[:solid_inlet_loading, Unit()])
end

"""
    SolidInletTemperature

Temperature of the sorbent entering the top of a unit that moves its sorbent
[K]. Not used when the top is fed by another unit.
"""
struct SolidInletTemperature <: Jutul.ScalarVariable end
Jutul.associated_entity(::SolidInletTemperature) = Unit()
function Jutul.default_parameter_values(data_domain, model, ::SolidInletTemperature, symb)
    return copy(data_domain[:solid_inlet_temperature, Unit()])
end

"""
    ParticlePorosity

Pore volume fraction of the sorbent particles of a unit that moves its
sorbent. The gas in the pores moves with the sorbent (see
[`PoreGasFlow`](@ref)); 0 carries none.
"""
struct ParticlePorosity <: Jutul.ScalarVariable end
Jutul.associated_entity(::ParticlePorosity) = Unit()
Jutul.minimum_value(::ParticlePorosity) = 0.0
function Jutul.default_parameter_values(data_domain, model, ::ParticlePorosity, symb)
    return copy(data_domain[:particle_porosity, Unit()])
end

"""
    SolidInletGasConcentration

Concentration of each component in the pores of the sorbent entering the top
of a unit that moves its sorbent [mol/m³]. Not used when the top is fed by
another unit.
"""
struct SolidInletGasConcentration <: ComponentVariable end
Jutul.associated_entity(::SolidInletGasConcentration) = Unit()
Jutul.minimum_value(::SolidInletGasConcentration) = 0.0
function Jutul.default_parameter_values(data_domain, model, ::SolidInletGasConcentration, symb)
    return copy(data_domain[:solid_inlet_gas_concentration, Unit()])
end

"Volume flow of sorbent particles [m³/s]."
@inline sorbent_volume_flow(state) = state.SorbentMassFlow[1] / state.AdsorbentDensity[1]

# Value of φ (a function of the cell) at the bottom face of `cell`, where the
# sorbent leaves it
@inline _sorbent_face(limiter, φ, cell, nc) = _face_value(limiter, φ, cell, cell < nc ? cell + 1 : nothing, cell > 1 ? cell - 1 : nothing)

# Values of φ where the sorbent enters and leaves `cell`; φ_inlet is the value
# entering the top cell
@inline function _sorbent_faces(limiter, φ, φ_inlet, cell, nc)
    φ_out = _sorbent_face(limiter, φ, cell, nc)
    φ_in = cell == nc ? φ_inlet() : _sorbent_face(limiter, φ, cell + 1, nc)
    return (φ_in, φ_out)
end

# Q_s·(q_out − q_in) in the loading equation of component i, if the sorbent moves
@inline function _loading_advection(model, state, i, cell)
    moves_sorbent(model.system) || return false
    nc = Jutul.number_of_cells(model.domain)
    q = state.AdsorbedLoading
    q_in, q_out = _sorbent_faces(advection_limiter(model.system), c -> q[i, c], () -> state.SolidInletLoading[i, 1], cell, nc)
    return sorbent_volume_flow(state) * (q_out - q_in)
end

"""
Gas taken up by moving sorbent, at the mass transfer rate: removes
`V_s·rate_i` from the gas mass balance. Of the gas adsorbed, `V_s·∂q/∂t`
stays in the cell and `Q_s·(q − q_above)` leaves with the sorbent.
"""
struct SorptionRate <: AbstractSourceTerm end

@inline function Framework.mass_source(::SorptionRate, model, state, state0, cell, i, Δt)
    return -state.SolidVolume[cell] * state.AdsorptionMassTransfer[i, cell]
end

"""
Heat released by sorption on moving sorbent, and the sensible heat that the
moles taken up carry from the gas to the adsorbed phase:
[`HeatOfSorption`](@ref) with the mass transfer rate in place of `∂q/∂t`.
Removes `V_s·Σ_i (C_pg·T + Δh_i + ΔH_i)·rate_i` from the gas energy balance.
"""
struct SorptionRateHeat <: AbstractSourceTerm end

@inline function Framework.energy_source(::SorptionRateHeat, model, state, state0, cell, Δt)
    rate = state.AdsorptionMassTransfer
    ΔH = state.ΔH
    T = state.Temperature[cell]
    gas = state.C_pg[cell] * T
    Δh = _adsorbed_excess_enthalpies(model.system, T)
    term = zero(T)
    for i in axes(rate, 1)
        term += (gas + Δh[i] + ΔH[i, cell]) * rate[i, cell]
    end
    return -state.SolidVolume[cell] * term
end

"""
Sensible heat brought into a cell by the sorbent moving through it:
`Q_s·(ρ_s·C_ps + C_pa·Σ_i q_i)·(T_in − T_out)`, with the sorbent
temperatures where it enters and leaves the cell, the moving-sorbent
counterpart of [`SolidHeatCapacity`](@ref) and
[`AdsorbedPhaseHeatCapacity`](@ref).
"""
struct SorbentHeatFlow <: AbstractSourceTerm end

"Heat capacity flow [W/K] of the sorbent leaving `cell`."
@inline function sorbent_heat_capacity_flow(state, cell)
    AC = state.AdsorbedConcentration
    sq = zero(eltype(AC))
    for i in axes(AC, 1)
        sq += AC[i, cell]
    end
    c_solid = state.AdsorbentDensity[1] * state.AdsorbentHeatCapacity[1]
    c_adsorbed = state.C_pa[cell] * sq
    return sorbent_volume_flow(state) * (c_solid + c_adsorbed)
end

@inline function Framework.energy_source(::SorbentHeatFlow, model, state, state0, cell, Δt)
    nc = Jutul.number_of_cells(model.domain)
    T = state.Temperature
    T_in, T_out = _sorbent_faces(advection_limiter(model.system), c -> T[c], () -> state.SolidInletTemperature[1], cell, nc)
    return sorbent_heat_capacity_flow(state, cell) * (T_in - T_out)
end

"""
Gas in the pores of the sorbent, carried down with it: the sorbent entering a
cell brings `Q_s·ε_p·c_in` of each component and its energy, and the
sorbent leaving takes `Q_s·ε_p·c_out`, with the concentrations and
temperatures where it enters and leaves. The pore gas has the composition and
temperature of the gas around the particles, and counts in the gas holdup of
the cell (`FluidVolume`). Its energy is in the form of the gas energy balance
(see [`pore_gas_energy`](@ref)).
"""
struct PoreGasFlow <: AbstractSourceTerm end

# Volume flow of pore gas [m³/s]
@inline pore_gas_volume_flow(state) = sorbent_volume_flow(state) * state.ParticlePorosity[1]


"""
    pore_gas_energy(sys, state, cell, c, T_gas) → e

Energy [J/m³ of pore gas] that pore gas with concentrations `c` [mol/m³] at
temperature `T_gas` [K] brings into `cell` of a unit. Per mole it is the
`C_pg·T` of the cell, which the gas energy terms count for every mole the
cell gains or loses (energy.jl), and the enthalpy that brings the gas from
the cell's temperature to `T_gas`: `Σ_i c_i·(C_pg·T + h_i(T_gas) − h_i(T))`,
as for a stream from a device (networks/ports.jl). Gas leaving at the
temperature of the cell takes `C_pg·T` per mole.
"""
@inline function pore_gas_energy(sys, state, cell, c, T_gas)
    T = state.Temperature[cell]
    Δh = component_enthalpies(sys.thermodynamics, T_gas) .- component_enthalpies(sys.thermodynamics, T)
    cT = state.C_pg[cell] * T
    return sum(i -> c[i] * (cT + Δh[i]), eachindex(Δh))
end

@inline function Framework.mass_source(::PoreGasFlow, model, state, state0, cell, i, Δt)
    nc = Jutul.number_of_cells(model.domain)
    c = state.MolarConcentration
    c_in, c_out = _sorbent_faces(advection_limiter(model.system), k -> c[i, k], () -> state.SolidInletGasConcentration[i, 1], cell, nc)
    return pore_gas_volume_flow(state) * (c_in - c_out)
end

@inline function Framework.energy_source(::PoreGasFlow, model, state, state0, cell, Δt)
    sys = model.system
    nc = Jutul.number_of_cells(model.domain)
    lim = advection_limiter(sys)
    c, T = state.MolarConcentration, state.Temperature
    faces = ntuple(i -> _sorbent_faces(lim, k -> c[i, k], () -> state.SolidInletGasConcentration[i, 1], cell, nc), number_of_components(sys))
    T_in, T_out = _sorbent_faces(lim, k -> T[k], () -> state.SolidInletTemperature[1], cell, nc)
    e_in = pore_gas_energy(sys, state, cell, map(first, faces), T_in)
    e_out = pore_gas_energy(sys, state, cell, map(last, faces), T_out)
    return pore_gas_volume_flow(state) * (e_in - e_out)
end
