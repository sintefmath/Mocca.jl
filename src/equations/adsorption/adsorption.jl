
function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:AdsorbedConcentration},
    model::AdsorptionModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    ϵ = model.system.p.Φ
    solid_volume = state[:SolidVolume]

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] =
            solid_volume[self_cell] * ∂M∂t -
            solid_volume[self_cell] * forcing_term[component, self_cell] 
    end
end
