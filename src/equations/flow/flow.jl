@inline function Jutul.face_flux!(
    q,
    face,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    state,
    model::AdsorptionModel,
    dt,
    flow_disc::Jutul.PotentialFlow,
    ldisc,
)
    kgrad, upw = ldisc.face_disc(face)

    c = state.MolarConcentration
    μ = state.FluidViscosity[1]

    T_f = state.Transmissibilities[face]
    ∇p = Jutul.gradient(state.Pressure, kgrad)
    q_darcy = -T_f * ∇p
    L = kgrad.left
    R = kgrad.right

    cL = state.TotalMolarConcentration[L]
    cR = state.TotalMolarConcentration[R]
    y = state.y
    C = (cL + cR)/2.0

    D_l = state.DiffusionTransmissibilities[face]
    for component in eachindex(q)
        F_c = cell -> c[component, cell] / μ
        c_face = Jutul.upwind(upw, F_c, q_darcy)
        q_i = c_face * q_darcy - C * D_l * Jutul.gradient(y, component, kgrad)

        q = setindex(q, q_i, component)
    end
    return q
end

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    model::AdsorptionModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    # Compute ∇⋅V
    disc = eq.flow_discretization
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)

    AC = state.AdsorbedConcentration
    AC₀ = state0.AdsorbedConcentration

    AC = state.AdsorbedConcentration
    AC₀ = state0.AdsorbedConcentration
    @inbounds V = state.SolidVolume[self_cell]
    @inbounds for i in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
        ∂q∂t = Jutul.accumulation_term(AC, AC₀, Δt, i, self_cell)
        V = state.SolidVolume[self_cell]
        eq_buf[i] = ∂M∂t + div_v[i] + (V * ∂q∂t)
    end
end

