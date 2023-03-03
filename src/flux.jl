function JutulDarcy.component_mass_fluxes!(
    q,
    face,
    state,
    model::Jutul.SimulationModel{G,S},
    kgrad,
    upw,
) where {G<:Any,S<:AdsorptionFlowSystem}
    # This is defined for us:
    # kgrad = TPFA(left, right, face_sign)
    # upw = SPU(left, right)
    # TODO: Implement me

    sys = model.system
    disc = JutulDarcy.kgrad_common(face, state, model, kgrad)
    (∇p, T_f, gΔz) = disc
    # @show T_f
    # @show ∇p
    c = state.concentrations
    μ = sys.fluid_viscosity
    q_darcy = -T_f * ∇p

    R = sys.R
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
    eq::Jutul.ConservationLaw{:TotalMasses},
    model::AdsorptionFlowModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    #@info "Updating total masses"

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    # Compute ∇⋅V
    disc = eq.flow_discretization
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)
    # @info "Forcing term" forcing_term
    forcing_term_coefficient = model.system.forcing_term_coefficient
    for i in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
        # @info i ∂M∂t forcing_term[i, self_cell]
        eq_buf[i] = ∂M∂t + div_v[i] + forcing_term_coefficient * forcing_term[i, self_cell]
    end
end


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
    ϵ = model.system.Φ
    # @info "From transfer system" forcing_term


    # TODO: Figure out a better way to compute solid volume
    solid_volume = state[:solidVolume]
    forcing_term_coefficient = model.system.forcing_term_coefficient

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] =
            ∂M∂t -
            forcing_term_coefficient * forcing_term[component, self_cell] /
            solid_volume[self_cell]
    end
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

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = ∂M∂t
    end
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

    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = ∂M∂t
    end
end