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

"""
Viscosity of the gas in a cell relative to the `FluidViscosity` parameter, from
the system's `transport` (a TransportReference). A unit whose transport
properties do not change leaves it out, and its viscosity is the parameter.
"""
struct ViscosityFactor <: Jutul.ScalarVariable end

Jutul.@jutul_secondary function update_viscosity_factor!(f, tv::ViscosityFactor, model::MoccaModel, y, Temperature, ix)
    transport = model.system.transport
    for cell in ix
        f[cell] = viscosity_factor(transport, view(y, :, cell), Temperature[cell])
    end
end

# Viscosity [Pa·s] of the gas in a cell
@inline _cell_viscosity(state, cell) =
    _state_has(state, :ViscosityFactor) ? state.FluidViscosity[1] * state.ViscosityFactor[cell] : state.FluidViscosity[1]

"""
Molecular diffusivity [m²/s] of each component in a cell, for a dispersion
that follows it (see [`WakaoFunazkri`](@ref)): the diffusivities of the
system's `dispersion` model, at the feed, scaled by its `transport`.
"""
struct MolecularDiffusivity <: ComponentVariable end

Jutul.@jutul_secondary function update_molecular_diffusivity!(D, tv::MolecularDiffusivity, model::MoccaModel, y, Temperature, Pressure, ix)
    sys = model.system
    for cell in ix
        d = _scaled_diffusivities(sys.transport, sys.dispersion.diffusivities, view(y, :, cell), Temperature[cell], Pressure[cell])
        for i in axes(D, 1)
            D[i, cell] = d[i]
        end
    end
end

_scaled_diffusivities(r::TransportReference, D, y, T, P) = D .* diffusivity_factors(r, y, T, P)
_scaled_diffusivities(::TransportReference{ConstantTransport}, D, y, T, P) = D .+ zero(T / P)

"""
    DispersionScale

Factor on the axial dispersion of a unit whose dispersion follows the state
(see [`WakaoFunazkri`](@ref)), 1 by default. As a parameter it can be set per
case, by a design variable or in a history match, as `AxialDispersion` can for
a constant dispersion. Such a unit has it in place of `AxialDispersion`.
"""
struct DispersionScale <: Jutul.ScalarVariable end
Jutul.associated_entity(::DispersionScale) = Unit()
Jutul.minimum_value(::DispersionScale) = 0.0
Jutul.default_parameter_values(data_domain, model, ::DispersionScale, symb) = [1.0]

# Dispersion conductance [m³/s] of a face for a component: the AxialDispersion
# parameter, or by Wakao and Funazkri from the molecular diffusivity of the
# component in the cells of the face and the superficial velocity `u` [m/s]
@inline function _face_dispersion(state, model, face, kgrad, u, component)
    _state_has(state, :MolecularDiffusivity) || return state.AxialDispersion[1] * state.DispersionTransmissibilities[face]
    D = state.MolecularDiffusivity
    D̄ = (D[component, kgrad.left] + D[component, kgrad.right]) / 2
    return state.DispersionScale[1] * dispersion_coefficient(model.system.dispersion, D̄, u) * state.DispersionTransmissibilities[face]
end

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

# Darcy advection with upwinded concentration over viscosity plus Fickian
# dispersion of the mole fractions, constant or by Wakao and Funazkri.
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

    T_f = state.Transmissibilities[face]
    ∇p = Jutul.gradient(state.Pressure, kgrad)
    q_darcy = -T_f * ∇p
    L = kgrad.left
    R = kgrad.right

    cL = state.TotalMolarConcentration[L]
    cR = state.TotalMolarConcentration[R]
    y = state.y
    C = (cL + cR)/2.0

    # Superficial velocity, for a dispersion that follows it
    u = _state_has(state, :MolecularDiffusivity) ?
        q_darcy / (Jutul.upwind(upw, cell -> _cell_viscosity(state, cell), q_darcy) * state.BedCrossSectionArea[1]) : zero(q_darcy)
    for component in eachindex(q)
        F_c = cell -> c[component, cell] / _cell_viscosity(state, cell)
        c_face = Jutul.upwind(upw, F_c, q_darcy)
        D_c = _face_dispersion(state, model, face, kgrad, u, component)
        q_i = c_face * q_darcy - C * D_c * Jutul.gradient(y, component, kgrad)

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
