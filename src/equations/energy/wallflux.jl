
@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:WallConservedEnergy,<:Any},
    state,
    model::AdsorptionModel,
    dt,
    disc,
    flow_disc, 
    T = Float64
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)
    K_w = state.WallConductivity[1]
    T = state.WallTemperature
    q = K_w * JutulDarcy.gradient(T, kgrad)
    return q
end

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    model::AdsorptionModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    disc = eq.flow_discretization
    flux_temp(face) = face_flux_temperature(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_temp = ldisc.div(flux_temp)

    T = state.Temperature[self_cell]
    T_w = state.WallTemperature[self_cell]
    aw_in = state.WallAreaIn[self_cell]
    aw_out = state.WallAreaOut[self_cell]
    Δx = state.CellDx[self_cell]

    h_in = state.InnerHeatTransferCoeff[1]
    h_out = state.OuterHeatTransferCoeff[1]
    T_a = state.AmbientTemperature[1]
    C_pw = state.WallHeatCapacity[1]
    ρ_w = state.WallDensity[1]

    # This is from the paper:
    #source_term = 2 * r_in*h_in / (r_out^2-r_in^2)*(T-T_w) - 2 * r_out*h_out/(r_out^2-r_in^2) * (T_w - T_a)
    #this is from the matlab code:
    source_term = aw_in * h_in * (T-T_w) - aw_out * h_out * (T_w - T_a)

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        A_w = state.WallCrossSectionArea[1]
        wall_volume = A_w * Δx
        eq_buf[component] = wall_volume * ρ_w * C_pw * ∂M∂t - A_w * div_temp / Δx - source_term
    end
end
