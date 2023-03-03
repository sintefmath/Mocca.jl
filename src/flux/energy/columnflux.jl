@inline function Jutul.face_flux!(
    q_i,
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::AdsorptionFlowModel,
    dt,
    disc,
    flow_disc, 
    T = Float64
)
    # Specific version for tpfa flux
    # kgrad = TPFA(left, right, face_sign)
    #upw = SPU(left, right)
    #return component_mass_fluxes!(q_i, face, state, model, kgrad, upw)
    #@info "In face_flux!" q_i
    return q_i
end


function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    model::AdsorptionFlowModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    disc = eq.flow_discretization
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = ∂M∂t + div_v[component]
    end
end