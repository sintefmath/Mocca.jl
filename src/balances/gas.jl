# Gas-phase block: secondary variables, the advective–dispersive flux and the
# per-component mass balance. Shared by every distributed unit with a gas phase.
#
# A system with a gas phase must have a `thermodynamics` field, an
# AbstractGasThermodynamics for its equation of state and heat capacity.

struct AverageMolarMass <: Jutul.ScalarVariable end

struct MolarConcentration <: ComponentVariable end

# Total molar concentration (scalar per cell, from the equation of state)
struct TotalMolarConcentration <: Jutul.ScalarVariable end
Jutul.minimum_value(::TotalMolarConcentration) = 0.0

# Per-component moles in cell (mol)
struct ComponentMasses <: ComponentVariable end
Jutul.minimum_value(::ComponentMasses) = 0.0

"Molar heat capacity [J/(mol·K)] of the gas in a cell, from the system's `thermodynamics`."
struct GasMolarHeatCapacity <: Jutul.ScalarVariable end

"""
    AxialDispersion

Axial dispersion coefficient [m²/s] of the gas in a distributed unit, one value
for the unit. It is a parameter, so it can be changed per case, set by a
design variable and differentiated without rebuilding the model.
"""
struct AxialDispersion <: Jutul.ScalarVariable end
Jutul.associated_entity(::AxialDispersion) = Unit()
Jutul.minimum_value(::AxialDispersion) = 0.0

function Jutul.default_parameter_values(data_domain, model, ::AxialDispersion, symb)
    haskey(data_domain, :diffusion_coefficient, Unit()) ||
        error(":diffusion_coefficient on Unit() must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    return copy(data_domain[:diffusion_coefficient, Unit()])
end

# Face transmissibilities [m] of the pore space for dispersion, per unit of
# AxialDispersion: the dispersive flux through a face is
# AxialDispersion·DispersionTransmissibilities·c̄·Δy
struct DispersionTransmissibilities <: Jutul.ScalarVariable end
Jutul.minimum_value(::DispersionTransmissibilities) = 0.0
Jutul.associated_entity(::DispersionTransmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::DispersionTransmissibilities, symb)
    g = Jutul.physical_representation(data_domain)
    return Jutul.compute_face_trans(g, data_domain[:porosity])
end

Jutul.@jutul_secondary function update_total_concentration!(ctot, tv::TotalMolarConcentration, model::MoccaModel, Pressure, Temperature, ix)
    thermo = model.system.thermodynamics
    for cell in ix
        ctot[cell] = molar_concentration(thermo, Pressure[cell], Temperature[cell])
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

Jutul.@jutul_secondary function update_heat_capacity_fluid!(cpg, tv::GasMolarHeatCapacity, model::MoccaModel, y, Temperature, ix)
    thermo = model.system.thermodynamics
    for cell in ix
        cpg[cell] = molar_heat_capacity(thermo, view(y, :, cell), Temperature[cell])
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

    D_l = state.AxialDispersion[1] * state.DispersionTransmissibilities[face]
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
