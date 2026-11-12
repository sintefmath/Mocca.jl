# Heat exchanger block: tubes through a distributed unit carrying a heat
# transfer fluid, for units with the HeatExchanger thermal model. Per cell
# (Mondino et al. 2022, eqs. A.10–A.11, with areas and volumes per cell):
#
#   tube wall:  V_t ρ_t c_t ∂T_t/∂t = A_o h_o (T − T_t) − A_i h_i (T_t − T_f)
#   fluid:      V_f ρ_f c_f ∂T_f/∂t + ṁ_f c_f (T_f,out − T_f,in) = A_i h_i (T_t − T_f)
#
# and the gas–sorbent energy balance loses A_o h_o (T − T_t). The fluid flows
# from the bottom cell to the top cell; it enters the bottom cell at
# HXFluidInletTemperature, or from another unit through a
# HeatTransferFluidConnection (networks/sorbent.jl). The temperatures where it
# enters and leaves a cell are reconstructed with the system's limiter, first
# order at the faces above the bottom cell and the top cell (see
# solid_advection.jl).

struct TubeEnergy <: Energy end
struct HXFluidEnergy <: Energy end

# Secondary copies of the temperatures, for the time differences; see
# AdsorbedLoading in blocks/sorbent.jl for why.
Jutul.@jutul_secondary function update_tube_conserved_energy!(e, tv::TubeEnergy, model::DistributedUnitModel, TubeTemperature, ix)
    for cell in ix
        e[cell] = TubeTemperature[cell]
    end
end

Jutul.@jutul_secondary function update_hx_fluid_conserved_energy!(e, tv::HXFluidEnergy, model::DistributedUnitModel, HXFluidTemperature, ix)
    for cell in ix
        e[cell] = HXFluidTemperature[cell]
    end
end

# Per-cell parameters, from a per-volume value in the data domain
abstract type _PerVolumeParameter <: Jutul.ScalarVariable end
Jutul.minimum_value(::_PerVolumeParameter) = 0.0

function Jutul.default_parameter_values(data_domain, model, p::_PerVolumeParameter, symb)
    key = _domain_key(p)
    haskey(data_domain, key, Unit()) || error("$key on Unit() must be present in the DataDomain to initialize $symb; build the unit with a heat_exchanger")
    return first(data_domain[key, Unit()]) .* data_domain[:volumes]
end

"Outer tube area in each cell [m²], facing the gas and sorbent."
struct TubeOuterArea <: _PerVolumeParameter end
"Inner tube area in each cell [m²], facing the heat transfer fluid."
struct TubeInnerArea <: _PerVolumeParameter end
"Volume of tube wall material in each cell [m³]."
struct TubeVolume <: _PerVolumeParameter end
"Volume of heat transfer fluid in each cell [m³]."
struct HXFluidVolume <: _PerVolumeParameter end

_domain_key(::TubeOuterArea) = :hx_outer_area_density
_domain_key(::TubeInnerArea) = :hx_inner_area_density
_domain_key(::TubeVolume) = :hx_tube_volume_fraction
_domain_key(::HXFluidVolume) = :hx_fluid_volume_fraction

# Unit parameters, read from the data domain under the given key
abstract type _HXUnitParameter <: Jutul.ScalarVariable end
Jutul.associated_entity(::_HXUnitParameter) = Unit()

function Jutul.default_parameter_values(data_domain, model, p::_HXUnitParameter, symb)
    key = _domain_key(p)
    haskey(data_domain, key, Unit()) || error("$key on Unit() must be present in the DataDomain to initialize $symb; build the unit with a heat_exchanger")
    return copy(data_domain[key, Unit()])
end

"Density of the tube wall material [kg/m³]."
struct TubeDensity <: _HXUnitParameter end
"Specific heat capacity of the tube wall material [J/(kg·K)]."
struct TubeHeatCapacity <: _HXUnitParameter end
"Density of the heat transfer fluid [kg/m³]."
struct HXFluidDensity <: _HXUnitParameter end
"Specific heat capacity of the heat transfer fluid [J/(kg·K)]."
struct HXFluidHeatCapacity <: _HXUnitParameter end
"Heat transfer coefficient between the gas and sorbent and the outer tube surface [W/(m²·K)]."
struct SorbentSideHeatTransferCoeff <: _HXUnitParameter end
"Heat transfer coefficient between the inner tube surface and the fluid [W/(m²·K)]."
struct FluidSideHeatTransferCoeff <: _HXUnitParameter end
"Mass flow of heat transfer fluid through the unit [kg/s], from bottom to top."
struct HXFluidMassFlow <: _HXUnitParameter end
"Temperature of the heat transfer fluid entering the bottom of the unit [K]. Not used when it comes from another unit."
struct HXFluidInletTemperature <: _HXUnitParameter end

_domain_key(::TubeDensity) = :hx_tube_density
_domain_key(::TubeHeatCapacity) = :hx_tube_heat_capacity
_domain_key(::HXFluidDensity) = :hx_fluid_density
_domain_key(::HXFluidHeatCapacity) = :hx_fluid_heat_capacity
_domain_key(::SorbentSideHeatTransferCoeff) = :hx_sorbent_side_htc
_domain_key(::FluidSideHeatTransferCoeff) = :hx_fluid_side_htc
_domain_key(::HXFluidMassFlow) = :hx_fluid_mass_flow
_domain_key(::HXFluidInletTemperature) = :hx_fluid_inlet_temperature

Jutul.minimum_value(::HXFluidMassFlow) = 0.0

function select_heat_exchanger_parameters!(S)
    S[:TubeOuterArea] = TubeOuterArea()
    S[:TubeInnerArea] = TubeInnerArea()
    S[:TubeVolume] = TubeVolume()
    S[:HXFluidVolume] = HXFluidVolume()
    S[:TubeDensity] = TubeDensity()
    S[:TubeHeatCapacity] = TubeHeatCapacity()
    S[:HXFluidDensity] = HXFluidDensity()
    S[:HXFluidHeatCapacity] = HXFluidHeatCapacity()
    S[:SorbentSideHeatTransferCoeff] = SorbentSideHeatTransferCoeff()
    S[:FluidSideHeatTransferCoeff] = FluidSideHeatTransferCoeff()
    S[:HXFluidMassFlow] = HXFluidMassFlow()
    S[:HXFluidInletTemperature] = HXFluidInletTemperature()
    return S
end

# Heat [W] from the gas and sorbent to the tube wall, and from the wall to the fluid
@inline _outer_tube_heat(state, cell) = state.TubeOuterArea[cell] * state.SorbentSideHeatTransferCoeff[1] * (state.Temperature[cell] - state.TubeTemperature[cell])
@inline _inner_tube_heat(state, cell) = state.TubeInnerArea[cell] * state.FluidSideHeatTransferCoeff[1] * (state.TubeTemperature[cell] - state.HXFluidTemperature[cell])

"Heat capacity flow of the heat transfer fluid [W/K]."
@inline hx_fluid_heat_capacity_flow(state) = state.HXFluidMassFlow[1] * state.HXFluidHeatCapacity[1]

"Heat exchange between the gas and sorbent and the tubes of a heat exchanger, `A_o·h_o·(T − T_t)`."
struct TubeExchange <: AbstractSourceTerm end

@inline Framework.energy_source(::TubeExchange, model, state, state0, cell, Δt) = -_outer_tube_heat(state, cell)

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:TubeConservedEnergy},
    model::DistributedUnitModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    Δt = _scaled_dt(state, Δt)
    ∂T∂t = Jutul.accumulation_term(state.TubeConservedEnergy, state0.TubeConservedEnergy, Δt, 1, self_cell)
    C = state.TubeVolume[self_cell] * state.TubeDensity[1] * state.TubeHeatCapacity[1]
    eq_buf[1] = C * ∂T∂t - _outer_tube_heat(state, self_cell) + _inner_tube_heat(state, self_cell)
end

# Fluid temperature at the top face of `cell`, where the fluid leaves it
@inline _hx_fluid_face(limiter, T, cell, nc) = _face_value(limiter, c -> T[c], cell, cell > 1 ? cell - 1 : nothing, cell < nc ? cell + 1 : nothing)

function Jutul.update_equation_in_entity!(
    eq_buf::AbstractVector{T_e},
    self_cell,
    state,
    state0,
    eq::Jutul.ConservationLaw{:HXFluidConservedEnergy},
    model::DistributedUnitModel,
    Δt,
    ldisc = Jutul.local_discretization(eq, self_cell),
) where {T_e}
    Δt = _scaled_dt(state, Δt)
    ∂T∂t = Jutul.accumulation_term(state.HXFluidConservedEnergy, state0.HXFluidConservedEnergy, Δt, 1, self_cell)
    C = state.HXFluidVolume[self_cell] * state.HXFluidDensity[1] * state.HXFluidHeatCapacity[1]
    T = state.HXFluidTemperature
    nc = Jutul.number_of_cells(model.domain)
    limiter = advection_limiter(model.system)
    T_out = _hx_fluid_face(limiter, T, self_cell, nc)
    T_in = self_cell == 1 ? state.HXFluidInletTemperature[1] : _hx_fluid_face(limiter, T, self_cell - 1, nc)
    advection = hx_fluid_heat_capacity_flow(state) * (T_out - T_in)
    eq_buf[1] = C * ∂T∂t + advection - _inner_tube_heat(state, self_cell)
end

# Convergence, as for the other energy balances (convergence.jl): the residual
# as the temperature error it leaves over a step. The heat capacity of the
# tube wall or fluid in a cell may be small next to the heat it exchanges, so
# the measure divides by C/Δt plus the conductances to its neighbours.
function Jutul.convergence_criterion(model::DistributedUnitModel, storage, eq::Jutul.ConservationLaw{:TubeConservedEnergy}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    function per_kelvin(cell)
        C = state.TubeVolume[cell] * state.TubeDensity[1] * state.TubeHeatCapacity[1]
        UA = state.TubeOuterArea[cell] * state.SorbentSideHeatTransferCoeff[1] + state.TubeInnerArea[cell] * state.FluidSideHeatTransferCoeff[1]
        return Jutul.value(C + dt * UA)
    end
    e = _cell_criterion(per_kelvin, r, 1, dt)
    return (ΔT = (errors = e, names = ["ΔT_tube"]),)
end

function Jutul.convergence_criterion(model::DistributedUnitModel, storage, eq::Jutul.ConservationLaw{:HXFluidConservedEnergy}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    function per_kelvin(cell)
        C = state.HXFluidVolume[cell] * state.HXFluidDensity[1] * state.HXFluidHeatCapacity[1]
        UA = state.TubeInnerArea[cell] * state.FluidSideHeatTransferCoeff[1] + hx_fluid_heat_capacity_flow(state)
        return Jutul.value(C + dt * UA)
    end
    e = _cell_criterion(per_kelvin, r, 1, dt)
    return (ΔT = (errors = e, names = ["ΔT_fluid"]),)
end
