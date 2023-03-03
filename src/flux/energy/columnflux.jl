@inline function Jutul.face_flux!(
    q,
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
    kgrad, upw = flow_disc.face_disc(face)

    sys = model.system
    disc = JutulDarcy.kgrad_common(face, state, model, kgrad)
    (∇p, T_f, gΔz) = disc
    # @show T_f
    # @show ∇p
    c = state.concentrations
    μ = sys.p.fluid_viscosity
    q_darcy = -T_f * ∇p

    R = sys.p.R
    L = kgrad.left
    R = kgrad.right

    favg(X) = (X[L] + X[R]) / 2
    P = favg(state.Pressure)
    T = favg(state.Temperature)
    C = P / (R * T)

    D_l = axial_dispersion(sys)
    for component in eachindex(q)
        F_c = cell -> c[component, cell] / μ
        c_face = JutulDarcy.upwind(upw, F_c, q_darcy)
        y_i = view(state.y, component, :)
        q_i = c_face * q_darcy + C * D_l * JutulDarcy.gradient(y_i, kgrad)

        q = setindex(q, q_i, component)
    end
    # @info "Flux $face" q
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
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = ∂M∂t + div_v[component]
    end
end