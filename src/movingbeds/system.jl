"""
    MovingBed{N, RealT, IsoT, MtT, ThT, LimT} <: SorbentBed

A moving bed of sorbent: the sorbent moves down through the unit at a steady
rate, from the top cell to the bottom cell, while gas flows through it in
either direction. It is the building block of moving bed temperature swing
adsorption (MBTSA), where the sorbent circulates through sections at different
temperatures (see [`SorbentConnection`](@ref) and [`add_moving_bed!`](@ref)).

Compared with a [`FixedBed`](@ref), the loading and the sorbent temperature are
carried with the sorbent: the loading equation gains the upwind term
`Q_s·(q − q_above)` and the energy balance the sensible heat of the moving
sorbent. Here `Q_s` is the volume flow of sorbent particles, the
`SorbentMassFlow` parameter divided by the particle density. The sorbent
entering the top cell has the `SolidInletLoading` and `SolidInletTemperature`
parameters, unless the top is fed by another unit through a
[`SorbentConnection`](@ref).

Sorption enters the gas balances through the mass transfer rate itself, rather
than through `∂q/∂t` as in a fixed bed, since part of the change in loading is
carried away with the sorbent. With no sorbent flow the two forms agree once
the loading equation is solved.

Gas and sorbent share one temperature, as in the fixed bed. The gas in the
particles' pores (a fraction `ParticlePorosity` of the particle volume) moves
with the sorbent and counts in the gas held in each cell; see
[`PoreGasFlow`](@ref).

Thermal models: [`Adiabatic`](@ref) (default), [`Isothermal`](@ref),
[`WithWall`](@ref) and [`HeatExchanger`](@ref).

What the sorbent carries, and the heat transfer fluid of a heat exchanger, is
advected with a limited second-order upwind scheme, [`VanAlbada`](@ref) by
default (see [`AbstractLimiter`](@ref)); [`Upwind`](@ref) is first order.
"""
struct MovingBed{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer, ThT<:AbstractThermalModel, LimT<:AbstractLimiter,
        GasT<:AbstractGasThermodynamics, AdsT<:AbstractGasThermodynamics} <: SorbentBed
    component_names::Vector{String}
    molecular_masses::SVector{N, RealT}         # Molecular masses per component [kg/mol]
    thermodynamics::GasT                        # Heat capacities and equation of state of the gas
    adsorbed_thermodynamics::AdsT               # Heat capacities of the adsorbed phase
    isotherm::IsoT
    mass_transfer::MtT
    thermal::ThT
    limiter::LimT
end

const MovingBedModel = Jutul.SimulationModel{<:Any, <:MovingBed, <:Any, <:Any}

"""
    MovingBed(; isotherm, mass_transfer, molecular_masses, component_names,
        thermodynamics, adsorbed_thermodynamics = thermodynamics,
        thermal = Adiabatic(), limiter = VanAlbada())

Construct an N-component moving bed system from explicit physics objects.
`thermodynamics` and `adsorbed_thermodynamics` are
[`AbstractGasThermodynamics`](@ref) models of the gas and the adsorbed phase;
instead of them, `heat_capacity_gas` and `heat_capacity_adsorbed` give
constant heat capacities per mass [J/(kg·K)], as for a
[`FixedBed`](@ref). The viscosity and diffusivities of the gas are those at
the feed. `limiter` reconstructs what the sorbent and the heat transfer fluid carry at
the cell faces (see [`AbstractLimiter`](@ref)).
"""
function MovingBed(;
    isotherm::AbstractIsotherm,
    mass_transfer::AbstractMassTransfer,
    molecular_masses,
    component_names,
    thermodynamics = nothing,
    adsorbed_thermodynamics = nothing,
    heat_capacity_gas = nothing,
    heat_capacity_adsorbed = nothing,
    thermal::AbstractThermalModel = Adiabatic(),
    limiter::AbstractLimiter = VanAlbada(),
)
    N = length(component_names)
    @assert length(molecular_masses) == N "molecular_masses length ($(length(molecular_masses))) must match component_names length ($N)"
    gas = _thermodynamics(thermodynamics, heat_capacity_gas, molecular_masses, "gas")
    isnothing(gas) && error("MovingBed needs the thermodynamics of the gas, or its heat_capacity_gas")
    adsorbed = something(_thermodynamics(adsorbed_thermodynamics, heat_capacity_adsorbed, molecular_masses, "adsorbed"), gas)
    @assert number_of_components(gas) == N "the gas thermodynamics must have $N components"
    @assert number_of_components(adsorbed) == N "the adsorbed thermodynamics must have $N components"
    RealT = eltype(molecular_masses)
    return MovingBed{N, RealT, typeof(isotherm), typeof(mass_transfer), typeof(thermal), typeof(limiter), typeof(gas), typeof(adsorbed)}(
        collect(String, component_names),
        SVector{N, RealT}(molecular_masses),
        gas,
        adsorbed,
        isotherm,
        mass_transfer,
        thermal,
        limiter,
    )
end

"""
    MovingBed(sys::MovingBed; thermal, limiter)

A copy of `sys` with another thermal model or limiter, for the sections of
one process that share a sorbent but are heated differently.
"""
function MovingBed(sys::MovingBed; thermal::AbstractThermalModel = sys.thermal, limiter::AbstractLimiter = sys.limiter)
    return MovingBed(;
        isotherm = sys.isotherm,
        mass_transfer = sys.mass_transfer,
        molecular_masses = sys.molecular_masses,
        component_names = sys.component_names,
        thermodynamics = sys.thermodynamics,
        adsorbed_thermodynamics = sys.adsorbed_thermodynamics,
        thermal = thermal,
        limiter = limiter)
end

# Terms the shared conservation laws pick up for a moving bed
Framework.moves_sorbent(::MovingBed) = true
Balances.advection_limiter(sys::MovingBed) = sys.limiter
Framework.mass_source_terms(::MovingBed) = (SorptionRate(), PoreGasFlow())
Framework.heat_capacity_terms(::MovingBed) = (SolidHeatCapacity(), AdsorbedPhaseHeatCapacity())
function Framework.energy_source_terms(sys::MovingBed)
    base = (SorptionRateHeat(), SorbentHeatFlow(), PoreGasFlow())
    if has_wall(sys)
        return (base..., WallExchange())
    elseif has_heat_exchanger(sys)
        return (base..., TubeExchange())
    else
        return base
    end
end
Framework.wall_source_terms(::MovingBed) = (WallExchange(), AmbientLoss(), WallEndConduction())

function Jutul.select_primary_variables!(S, system::MovingBed, model::Jutul.SimulationModel)
    S[:Pressure] = Pressure()
    S[:y] = GasMoleFractions()
    S[:AdsorbedConcentration] = AdsorbedConcentration()
    S[:Temperature] = Temperature()
    if has_wall(system)
        S[:WallTemperature] = Temperature()
    end
    if has_heat_exchanger(system)
        S[:TubeTemperature] = Temperature()
        S[:HXFluidTemperature] = Temperature()
    end
end

function Jutul.select_secondary_variables!(S, system::MovingBed, model::Jutul.SimulationModel)
    S[:TotalMolarConcentration] = TotalMolarConcentration()
    S[:MolarConcentration] = MolarConcentration()
    S[:AverageMolarMass] = AverageMolarMass()
    S[:ComponentMasses] = ComponentMasses()
    S[:AdsorptionMassTransfer] = AdsorptionMassTransfer()
    S[:AdsorbedLoading] = AdsorbedLoading()

    S[:ColumnConservedEnergy] = ColumnEnergy()
    if has_energy_balance(system)
        S[:ConservedPressure] = ConservedPressure()
    end
    if has_wall(system)
        S[:WallConservedEnergy] = WallEnergy()
    end
    if has_heat_exchanger(system)
        S[:TubeConservedEnergy] = TubeEnergy()
        S[:HXFluidConservedEnergy] = HXFluidEnergy()
    end
    S[:ΔH] = EnthalpyChange()
    S[:C_pa] = AdsorbedMolarHeatCapacity()
    S[:C_pg] = GasMolarHeatCapacity()
end

function Jutul.select_equations!(eqs, sys::MovingBed, model::Jutul.SimulationModel)
    fdisc = model.domain.discretizations.mass_flow
    nc = number_of_components(sys)

    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :ComponentMasses, nc)
    eqs[:mass_transfer] = Jutul.ConservationLaw(fdisc, :AdsorbedLoading, nc)
    if has_energy_balance(sys)
        eqs[:energy_column] = Jutul.ConservationLaw(fdisc, :ColumnConservedEnergy, 1)
    else
        eqs[:energy_column] = FixedTemperatureEquation()
    end
    if has_wall(sys)
        eqs[:energy_wall] = Jutul.ConservationLaw(fdisc, :WallConservedEnergy, 1)
    end
    if has_heat_exchanger(sys)
        eqs[:energy_tube] = Jutul.ConservationLaw(fdisc, :TubeConservedEnergy, 1)
        eqs[:energy_hx_fluid] = Jutul.ConservationLaw(fdisc, :HXFluidConservedEnergy, 1)
    end
end

function Jutul.select_parameters!(S, sys::MovingBed, model::Jutul.SimulationModel)
    # Per-cell parameters
    S[:SolidVolume] = BulkVolume()
    S[:FluidVolume] = FluidVolume()
    S[:CellDx] = CellDx()
    S[:Permeability] = Permeability()

    # Face transmissibilities
    S[:ThermalConductivities] = ThermalConductivities()
    S[:DispersionTransmissibilities] = DispersionTransmissibilities()

    # Unit-entity parameters
    S[:AxialDispersion] = AxialDispersion()
    S[:AdsorbentDensity] = AdsorbentDensity()
    S[:AdsorbentHeatCapacity] = AdsorbentHeatCapacity()
    S[:FluidViscosity] = FluidViscosity()
    S[:BedCrossSectionArea] = BedCrossSectionArea()
    S[:TimeScale] = TimeScale()

    # The moving sorbent
    S[:SorbentMassFlow] = SorbentMassFlow()
    S[:SolidInletLoading] = SolidInletLoading()
    S[:SolidInletTemperature] = SolidInletTemperature()
    S[:ParticlePorosity] = ParticlePorosity()
    S[:SolidInletGasConcentration] = SolidInletGasConcentration()

    if has_wall(sys)
        S[:WallAreaOut] = WallArea{:out}()
        S[:WallAreaIn] = WallArea{:in}()
        S[:WallDensity] = WallDensity()
        S[:WallHeatCapacity] = WallHeatCapacity()
        S[:WallConductivity] = WallConductivity()
        S[:InnerHeatTransferCoeff] = InnerHeatTransferCoeff()
        S[:OuterHeatTransferCoeff] = OuterHeatTransferCoeff()
        S[:AmbientTemperature] = AmbientTemperature()
        S[:WallCrossSectionArea] = WallCrossSectionArea()
    end
    if has_heat_exchanger(sys)
        select_heat_exchanger_parameters!(S)
    end
    if !has_energy_balance(sys)
        S[:IsothermalTemperature] = IsothermalTemperature()
    end
end
