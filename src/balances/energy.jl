# Gas-phase energy block, written in temperature form:
#
#   R = C_eff·∂T/∂t + (C_pg/R)·V_f·∂P/∂t − ∇·(K∇T) + (C_pg/R)·∇·(vP) − Σ energy sources
#
# where C_eff sums the unit's heat_capacity_terms and C_pg is the molar heat
# capacity of the gas in the cell. With P = cRT the gas terms are
# C_pg·(V_f·c·∂T/∂t + flow·∇T) plus C_pg·T for every mole the cell gains;
# the sources and streams that move moles in or out give that back (see
# HeatOfSorption and the stream cross term), so the balance holds for heat
# capacities that change with temperature and composition. Shared by every
# distributed unit with an energy balance.

abstract type Energy <: Jutul.ScalarVariable end
struct ColumnEnergy <: Energy end

# Pressure as a secondary variable, for the time difference in the pressure
# term; see AdsorbedLoading in balances/sorbent.jl for why.
struct ConservedPressure <: Jutul.ScalarVariable end

Jutul.@jutul_secondary function update_conserved_pressure!(p, tv::ConservedPressure, model::DistributedUnitModel, Pressure, ix)
    for cell in ix
        p[cell] = Pressure[cell]
    end
end

struct ThermalConductivities <: Jutul.ScalarVariable end
Jutul.variable_scale(::ThermalConductivities) = 1e-10
Jutul.minimum_value(::ThermalConductivities) = 0.0
Jutul.default_value(model, ::ThermalConductivities) = 1e-3
Jutul.associated_entity(::ThermalConductivities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalConductivities, symb)
    if haskey(data_domain, :thermal_conductivity, Unit())
        K = first(data_domain[:thermal_conductivity, Unit()])
        g = Jutul.physical_representation(data_domain)
        nc = Jutul.number_of_cells(g)
        T = Jutul.compute_face_trans(g, fill(K, nc))
    else
        error(":thermal_conductivity on Unit() must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

"""
Thermal conductivity of the gas in a cell relative to that at the reference of
the system's `transport`. The conductivity of the bed, the
`ThermalConductivities` of its faces, is scaled by it: the stagnant
conductivity of a packed bed is nearly proportional to that of its gas when
the solid conducts much better. A unit whose transport properties do not
change leaves it out.
"""
struct ConductivityFactor <: Jutul.ScalarVariable end

Jutul.@jutul_secondary function update_conductivity_factor!(f, tv::ConductivityFactor, model::DistributedUnitModel, y, Temperature, ix)
    transport = model.system.transport
    for cell in ix
        f[cell] = conductivity_factor(transport, view(y, :, cell), Temperature[cell])
    end
end

# Conductance of a face for heat: that of the parameter, scaled by the harmonic
# mean of the conductivity factors of its cells, as for conduction in series
@inline function _face_conductance(state, face, kgrad)
    K = state.ThermalConductivities[face]
    _state_has(state, :ConductivityFactor) || return K
    fL, fR = state.ConductivityFactor[kgrad.left], state.ConductivityFactor[kgrad.right]
    return K * (2 * fL * fR / (fL + fR))
end

Jutul.@jutul_secondary function update_column_conserved_energy!(column_energy, tv::ColumnEnergy, model::DistributedUnitModel, Temperature, ix)
    for cell in ix
        column_energy[cell] = Temperature[cell]
    end
end

@inline function face_flux_temperature(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::DistributedUnitModel,
    dt,
    disc,
    flow_disc,
    T=Float64,
)
    q = zero(Jutul.flux_vector_type(eq, T))
    kgrad, upw = flow_disc.face_disc(face)
    T = state.Temperature
    q = _face_conductance(state, face, kgrad) * Jutul.gradient(T, kgrad)
    return q
end

@inline function face_flux_pressure(
    face,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy,<:Any},
    state,
    model::DistributedUnitModel,
    dt,
    disc,
    flow_disc,
    T=Float64,
)
    q = zero(Jutul.flux_vector_type(eq, T))

    kgrad, upw = flow_disc.face_disc(face)

    T_f = state.Transmissibilities[face]
    ∇p = Jutul.gradient(state.Pressure, kgrad)
    # Upwinded as the concentration over viscosity of the gas flux
    μ = Jutul.upwind(upw, cell -> _cell_viscosity(state, cell), -T_f * ∇p)
    v = -T_f * ∇p / μ
    P_c = cell -> state.Pressure[cell]
    P_face = Jutul.upwind(upw, P_c, v)
    q = v * P_face
    return q
end

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    model::DistributedUnitModel,
    Δt,
    ldisc=Jutul.local_discretization(eq, self_cell),
) where {T_e}
    Δt = _scaled_dt(state, Δt)
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

    sys = model.system
    coeff_pressure = state.C_pg[self_cell] / GAS_CONSTANT

    ∂P∂t = (state.ConservedPressure[self_cell] - state0.ConservedPressure[self_cell]) / Δt
    pressure_term = coeff_pressure * state.FluidVolume[self_cell] * ∂P∂t

    accumulation_coeff = sum_terms(heat_capacity, heat_capacity_terms(sys), model, state, self_cell)
    src = sum_terms(energy_source, energy_source_terms(sys), model, state, state0, self_cell, Δt)

    for component in eachindex(eq_buf)
        ∂T∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        @inbounds eq_buf[component] =
            accumulation_coeff * ∂T∂t + pressure_term - div_temp +
            coeff_pressure * div_pressure[component] - src
    end
end
