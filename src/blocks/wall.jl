# Wall block: energy balance for a cylindrical wall around a distributed unit,
#
#   R = V_w·ρ_w·C_pw·∂T_w/∂t − A_w·∇·(K_w∇T_w)/Δx − Σ wall sources
#
# and the heat exchange terms that couple it to the gas and the surroundings.
# Used when the unit's thermal model is WithWall.

struct WallEnergy <: Energy end

Jutul.@jutul_secondary function update_wall_conserved_energy!(wall_energy, tv::WallEnergy, model::DistributedUnitModel, WallTemperature, ix)
    for cell in ix
        wall_energy[cell] = WallTemperature[cell]
    end
end

@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:WallConservedEnergy,<:Any},
    state,
    model::DistributedUnitModel,
    dt,
    disc,
    flow_disc,
    T = Float64
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)
    K_w = state.WallConductivity[1]
    T = state.WallTemperature
    q = K_w * Jutul.gradient(T, kgrad)
    return q
end

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    model::DistributedUnitModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    disc = eq.flow_discretization
    flux_temp(face) = face_flux_temperature(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_temp = ldisc.div(flux_temp)

    Δx = state.CellDx[self_cell]
    C_pw = state.WallHeatCapacity[1]
    ρ_w = state.WallDensity[1]
    A_w = state.WallCrossSectionArea[1]
    wall_volume = A_w * Δx

    src = sum_terms(wall_source, wall_source_terms(model.system), model, state, state0, self_cell, Δt)

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = wall_volume * ρ_w * C_pw * ∂M∂t - A_w * div_temp / Δx - src
    end
end

"""
Heat exchange between the gas and the inner wall surface, `A_in·h_in·(T − T_w)`.
It cools the gas and heats the wall.
"""
struct WallExchange <: AbstractSourceTerm end

@inline function _wall_exchange(state, cell)
    return state.WallAreaIn[cell] * state.InnerHeatTransferCoeff[1] * (state.Temperature[cell] - state.WallTemperature[cell])
end

@inline energy_source(::WallExchange, model, state, state0, cell, Δt) = -_wall_exchange(state, cell)
@inline wall_source(::WallExchange, model, state, state0, cell, Δt) = _wall_exchange(state, cell)

"Heat loss from the outer wall surface to ambient, `A_out·h_out·(T_w − T_a)`."
struct AmbientLoss <: AbstractSourceTerm end

@inline function wall_source(::AmbientLoss, model, state, state0, cell, Δt)
    return -state.WallAreaOut[cell] * state.OuterHeatTransferCoeff[1] * (state.WallTemperature[cell] - state.AmbientTemperature[1])
end

"""
Conduction along the wall out through its two ends to ambient, through a
half-cell connection at the first and last cell.
"""
struct WallEndConduction <: AbstractSourceTerm end

@inline function wall_source(::WallEndConduction, model, state, state0, cell, Δt)
    nc = Jutul.number_of_cells(model.domain)
    n_ends = (cell == 1) + (cell == nc)
    trans_wall = calc_bc_wall_trans(model, state, cell)
    return -n_ends * trans_wall * (state.WallTemperature[cell] - state.AmbientTemperature[1])
end

"""
    calc_bc_wall_trans(model, state, cell)

Thermal transmissibility [W/K] of the half-cell connection between wall `cell`
and the wall end it touches.
"""
function calc_bc_wall_trans(model, state, cell)
    k = state.WallConductivity[1]
    dx = state.CellDx[cell] / 2.0
    A = state.WallCrossSectionArea[1]
    return k * A / dx
end
