@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::AdsorptionModel,
    dt,
    disc,
    flow_disc,
    T=Float64,
)
    q = zero(Jutul.flux_vector_type(eq, T))
    kgrad, upw = flow_disc.face_disc(face)
    # K_z = model.system.p.K_z

    T = state.Temperature
    q = state.ThermalConductivities[face] * JutulDarcy.gradient(T, kgrad)
    return q
end

@inline function face_flux_pressure(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::AdsorptionModel,
    dt,
    disc,
    flow_disc,
    T=Float64,
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)

    sys = model.system
    T_f = JutulDarcy.effective_transmissibility(state, face, kgrad)
    ∇p = JutulDarcy.pressure_gradient(state, kgrad)
    μ = sys.p.fluid_viscosity
    v = -T_f * ∇p / μ
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
    model::AdsorptionModel,
    Δt,
    ldisc=Jutul.local_discretization(eq, self_cell),
) where {T_e}

    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]
    disc = eq.flow_discretization

    flux_temp(face) =
        face_flux_temperature(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    flux_pressure(face) =
        face_flux_pressure(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_temp = ldisc.div(flux_temp)
    div_pressure = ldisc.div(flux_pressure)

    C_pg = state.C_pg[self_cell]
    avm = state.avm[self_cell]

    T = state.Temperature[self_cell]
    T_w = state.WallTemperature[self_cell]

    A_win = state.WallAreaIn[self_cell]
    h_in = model.system.p.h_in
    R = model.system.p.R
    source_term = A_win * h_in * (T - T_w)
    ρ_s = model.system.p.ρ_s
    C_ps = model.system.p.C_ps
    C_pa = state.C_pa[self_cell]


    ∂P∂t = (state.Pressure[self_cell] - state0.Pressure[self_cell]) / Δt


    pv = state.FluidVolume
    sv = state.SolidVolume[self_cell]

    pressure_term = C_pg * avm / R * pv[self_cell] * ∂P∂t

    # HERE BE DRAGONS!
    #adsorption_term = state.SolidVolume[self_cell] * sum((state.C_pa[self_cell] * state.avm[self_cell] * state.Temperature[self_cell] .+ state.ΔH[:, self_cell]) .* ∂q∂t)
    # I think this should match now.
    # ∂q∂t =
    # (state.AdsorbedConcentration[:, self_cell] - state0.AdsorbedConcentration[:, self_cell]) ./ Δt

    # adsorption_term =
    #     sum((C_pa * avm * T * sv .* ∂q∂t)) + sum(state.ΔH[:, self_cell] * sv .* ∂q∂t)
    # sq = sum(state.AdsorbedConcentration[:, self_cell])
    # Rewritten as a loop:
    AC = state.AdsorbedConcentration
    AC₀ = state0.AdsorbedConcentration
    ΔH = state.ΔH

    sq = zero(T_e)
    adsorption_term = zero(T_e)
    for i in eachindex(eq_buf)
        ∂q∂t = Jutul.accumulation_term(AC, AC₀, Δt, i, self_cell)

        adsorption_term += (C_pa * avm * T + ΔH[i, self_cell]) * ∂q∂t
        sq += AC[i, self_cell]
    end
    # Finally multiply by volume
    adsorption_term *= sv

    accumulation_coeff = sv * (ρ_s * C_ps + C_pa * avm * sq)
    coeff_pressure = C_pg * avm / R

    for component in eachindex(eq_buf)
        #@info "Componennt" component size(eq_buf)
        ∂T∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)


        @inbounds eq_buf[component] =
            accumulation_coeff * ∂T∂t + pressure_term + adsorption_term - div_temp +
            coeff_pressure * div_pressure[component] +
            source_term[component]
    end
end
