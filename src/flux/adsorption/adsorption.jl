
function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:adsorptionRates},
    model::AdsorptionFlowModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    ϵ = model.system.p.Φ
    # @info "From transfer system" forcing_term


    # TODO: Figure out a better way to compute solid volume
    solid_volume = state[:solidVolume]
    forcing_term_coefficient = model.system.forcing_term_coefficient

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] =
            solid_volume[self_cell] * ∂M∂t -
            solid_volume[self_cell] * forcing_term_coefficient * forcing_term[component, self_cell] 
    end
end
