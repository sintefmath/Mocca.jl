# Gas-phase block: secondary variables, the advective–dispersive flux and the
# per-component mass balance. Shared by every distributed unit with a gas phase.

struct AverageMolarMass <: Jutul.ScalarVariable end

struct MolarConcentration <: ComponentVariable end

# Total molar concentration (scalar per cell, computed from ideal gas law)
struct TotalMolarConcentration <: Jutul.ScalarVariable end
Jutul.minimum_value(::TotalMolarConcentration) = 0.0

# Per-component moles in cell (mol)
struct ComponentMasses <: ComponentVariable end
Jutul.minimum_value(::ComponentMasses) = 0.0

struct SpecificHeatCapacityFluid <: Jutul.ScalarVariable end

struct DiffusionTransmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::DiffusionTransmissibilities) = 1e-10
Jutul.minimum_value(::DiffusionTransmissibilities) = 0.0
Jutul.default_value(model, ::DiffusionTransmissibilities) = 1e-3
Jutul.associated_entity(::DiffusionTransmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::DiffusionTransmissibilities, symb)
    if haskey(data_domain, :diffusion_coefficient, Unit())
        D = first(data_domain[:diffusion_coefficient, Unit()])
        g = Jutul.physical_representation(data_domain)
        porosity = data_domain[:porosity]
        T = Jutul.compute_face_trans(g, porosity .* D)
    else
        error(":diffusion_coefficient on Unit() must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

Jutul.@jutul_secondary function update_total_concentration!(ctot, tv::TotalMolarConcentration, model::MoccaModel, Pressure, Temperature, ix)
    for cell in ix
        ctot[cell] = Pressure[cell] / (GAS_CONSTANT * Temperature[cell])
    end
end

Jutul.@jutul_secondary function update_average_molar_mass!(avm, tv::AverageMolarMass, model::MoccaModel, y, ix)
    sys = model.system
    mm = sys.molecular_masses
    for cell in ix
        avm_val = zero(eltype(avm))
        for component in 1:number_of_components(sys)
            avm_val += y[component, cell] * mm[component]
        end
        avm[cell] = avm_val
    end
end

Jutul.@jutul_secondary function update_molar_concentration!(c, tv::MolarConcentration, model::MoccaModel, y, TotalMolarConcentration, ix)
    for cell in ix
        for component in 1:number_of_components(model.system)
            c[component, cell] = y[component, cell] * TotalMolarConcentration[cell]
        end
    end
end

Jutul.@jutul_secondary function update_component_masses!(m, tv::ComponentMasses, model::MoccaModel, MolarConcentration, FluidVolume, ix)
    for cell in ix
        for component in 1:number_of_components(model.system)
            m[component, cell] = MolarConcentration[component, cell] * FluidVolume[cell]
        end
    end
end

Jutul.@jutul_secondary function update_heat_capacity_fluid!(cpg, tv::SpecificHeatCapacityFluid, model::MoccaModel, y, ix)
    sys = model.system
    C_pg = sys.heat_capacity_gas
    T = eltype(cpg)
    for cell in ix
        cpg_i = zero(T)
        for component in 1:number_of_components(sys)
            cpg_i += y[component, cell] * C_pg[component]
        end
        cpg[cell] = cpg_i
    end
end

# Darcy advection with upwinded concentration plus Fickian dispersion of the
# mole fractions.
@inline function Jutul.face_flux!(
    q,
    face,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    state,
    model::DistributedUnitModel,
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

# R_i = ∂M_i/∂t + ∇·F_i − Σ mass sources
function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    model::DistributedUnitModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    Δt = _scaled_dt(state, Δt)
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    disc = eq.flow_discretization
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)

    terms = mass_source_terms(model.system)
    @inbounds for i in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
        src = sum_terms(mass_source, terms, model, state, state0, self_cell, i, Δt)
        eq_buf[i] = ∂M∂t + div_v[i] - src
    end
end
