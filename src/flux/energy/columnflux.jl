@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::AdsorptionFlowModel,
    dt,
    disc,
    flow_disc, 
    T = Float64
)
    q = zero(Jutul.flux_vector_type(eq, T))
    kgrad, upw = flow_disc.face_disc(face)
    K_z = model.system.p.K_z

    T = view(state.Temperature, :)
    q = K_z * JutulDarcy.gradient(T, kgrad)

    return q
end

@inline function face_flux_pressure(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::AdsorptionFlowModel,
    dt,
    disc,
    flow_disc, 
    T = Float64
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)

    sys = model.system
    disc = JutulDarcy.kgrad_common(face, state, model, kgrad)
    (∇p, T_f, gΔz) = disc
    v = -T_f * ∇p
    P_c = cell -> state.Pressure[cell]
    P_face = JutulDarcy.upwind(upw, P_c, v)
    q = v * P_face
    return q
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
    flux_temp(face) = face_flux_temperature(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    flux_pressure(face) = face_flux_pressure(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_temp = ldisc.div(flux_temp)
    div_pressure = ldisc.div(flux_pressure)

    C_pg = state.C_pg[self_cell]
    avm = state.avm[self_cell]

    T = state.Temperature[self_cell]
    T_w = state.WallTemperature[self_cell]
    A_w = area_wall(model.system)
    h_in = model.system.p.h_in
    R = model.system.p.R
    source_term = A_w * h_in * (T-T_w)
    ρ_s = model.system.p.ρ_s
    C_ps = model.system.p.C_ps
    C_pa = state.C_pa

    ∂P∂t = (state.Pressure[self_cell] - state0.Pressure[self_cell]) / Δt
    
    ∂q∂t = (state.adsorptionRates[:, self_cell] .- state0.adsorptionRates[:, self_cell]) ./ Δt

    pressure_term = C_pg * avm / R * ∂P∂t
    adsorption_term = state.solidVolume[self_cell] * sum((state.C_pa[self_cell] * state.avm[self_cell] * state.Temperature[self_cell] .+ state.ΔH[:, self_cell]) .* ∂q∂t)
    for component in eachindex(eq_buf)
        #@info "Componennt" component size(eq_buf)
        ∂T∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        sq = sum(state.adsorptionRates[:, self_cell])
        accumulation_coeff = state.solidVolume[self_cell] * (ρ_s * C_ps + C_pa[self_cell] * avm * sq)
        coeff_pressure = C_pg * avm / R
        eq_buf[component] = accumulation_coeff * ∂T∂t + pressure_term + adsorption_term - div_temp + coeff_pressure * div_pressure[component] + source_term[component]
    end
end