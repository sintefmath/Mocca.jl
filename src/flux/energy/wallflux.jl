
@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:WallConservedEnergy,<:Any},
    state,
    model::AdsorptionFlowModel,
    dt,
    disc,
    flow_disc, 
    T = Float64
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)
    K_w = model.system.p.K_w
    T = view(state.WallTemperature, :)
    q = K_w * JutulDarcy.gradient(T, kgrad)

    return q
end

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    model::AdsorptionFlowModel,
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
    h_in = model.system.p.h_in
    h_out = model.system.p.h_out
    r_in = model.system.p.r_in
    r_out = model.system.p.r_out
    T_a = model.system.p.T_a

    C_pw = model.system.p.C_pw
    ρ_w = model.system.p.ρ_w
    source_term = 2 * r_in*h_in / (r_out^2-r_in^2)*(T-T_w) - 2 * r_out*h_out/(r_out^2-r_in^2) * (T_w - T_a)

    
    Δx = compute_dx(model, self_cell)
    for component in eachindex(eq_buf)
        #@info "Componennt" component size(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        A_w = area_wall(model.system)
        wall_volume = A_w * Δx
        eq_buf[component] = wall_volume * ρ_w * C_pw * ∂M∂t - A_w * div_temp / Δx - wall_volume * source_term
    end
end