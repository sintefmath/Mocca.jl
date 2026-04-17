function Jutul.select_primary_variables!(
    S,
    system::AdsorptionSystem,
    model::Jutul.SimulationModel,
)
    S[:Pressure] = Pressure()
    S[:y] = GasMoleFractions()
    S[:AdsorbedConcentration] = AdsorbedConcentration()
    S[:Temperature] = Temperature()
    S[:WallTemperature] = Temperature()
end

function Jutul.select_secondary_variables!(
    S,
    system::AdsorptionSystem,
    model::Jutul.SimulationModel,
)
    S[:TotalMolarConcentration] = TotalMolarConcentration()
    S[:MolarConcentration] = MolarConcentration()
    S[:AverageMolarMass] = AverageMolarMass()
    S[:ComponentMasses] = ComponentMasses()
    S[:AdsorptionMassTransfer] = AdsorptionMassTransfer()


    # For the energy equations
    S[:ColumnConservedEnergy] = ColumnEnergy()
    S[:WallConservedEnergy] = WallEnergy()
    S[:ΔH] = EnthalpyChange()
    S[:C_pa] = SpecificHeatCapacityAdsorbent()
    S[:C_pg] = SpecificHeatCapacityFluid()
end

function Jutul.select_equations!(
    eqs,
    sys::AdsorptionSystem,
    model::Jutul.SimulationModel,
)
    fdisc = model.domain.discretizations.mass_flow
    nc = number_of_components(sys)

    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :ComponentMasses, nc)
    eqs[:mass_transfer] = Jutul.ConservationLaw(fdisc, :AdsorbedConcentration, nc)
    eqs[:energy_column] = Jutul.ConservationLaw(fdisc, :ColumnConservedEnergy, 1)
    eqs[:energy_wall] = Jutul.ConservationLaw(fdisc, :WallConservedEnergy, 1)
end

function Jutul.select_parameters!(S, ::AdsorptionSystem, model::Jutul.SimulationModel)
    # Per-cell parameters
    S[:SolidVolume] = BulkVolume()
    S[:FluidVolume] = FluidVolume()
    S[:WallAreaOut] = WallArea{:out}()
    S[:WallAreaIn] = WallArea{:in}()
    S[:CellDx] = CellDx()

    # Face transmissibilities
    S[:ThermalConductivities] = ThermalConductivities()
    S[:DiffusionTransmissibilities] = DiffusionTransmissibilities()

    # Column-entity parameters (single scalars)
    S[:AdsorbentDensity] = AdsorbentDensity()
    S[:AdsorbentHeatCapacity] = AdsorbentHeatCapacity()
    S[:WallDensity] = WallDensity()
    S[:WallHeatCapacity] = WallHeatCapacity()
    S[:WallConductivity] = WallConductivity()
    S[:FluidViscosity] = FluidViscosity()
    S[:FluidDensity] = FluidDensity()
    S[:InnerHeatTransferCoeff] = InnerHeatTransferCoeff()
    S[:OuterHeatTransferCoeff] = OuterHeatTransferCoeff()
    S[:AmbientTemperature] = AmbientTemperature()
    S[:WallCrossSectionArea] = WallCrossSectionArea()
end

