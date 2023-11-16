function Jutul.select_primary_variables!(
    S,
    system::AdsorptionSystem,
    model::Jutul.SimulationModel,
)
    S[:Pressure] = JutulDarcy.Pressure(minimum=π) # FIXME: Proper lower value 
    S[:y] = GasMoleFractions()
    S[:AdsorbedConcentration] = AdsorbedConcentration()
    S[:Temperature] = JutulDarcy.Temperature(min = 200.0, max_rel = 0.2)
    S[:WallTemperature] = JutulDarcy.Temperature(min = 200.0, max_rel = 0.2)
end

function Jutul.select_secondary_variables!(
    S,
    system::AdsorptionSystem,
    model::Jutul.SimulationModel,
)
    S[:cTot] = JutulDarcy.TotalMass()
    S[:concentrations] = Concentrations()
    S[:avm] = AverageMolecularMass()
    S[:TotalMasses] = JutulDarcy.TotalMasses()
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
    nc = JutulDarcy.number_of_components(sys)
    
    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :TotalMasses, nc)
    eqs[:mass_transfer] = Jutul.ConservationLaw(fdisc, :AdsorbedConcentration, nc)
    eqs[:energy_column] = Jutul.ConservationLaw(fdisc, :ColumnConservedEnergy, 1)
    eqs[:energy_wall] = Jutul.ConservationLaw(fdisc, :WallConservedEnergy, 1)
end

function Jutul.select_parameters!(S, ::AdsorptionSystem, model::Jutul.SimulationModel)
    S[:SolidVolume] = JutulDarcy.BulkVolume()
    S[:FluidVolume] = JutulDarcy.FluidVolume()
    S[:ThermalConductivities] = ThermalConductivities()
    S[:DiffusionTransmissibilities] = DiffusionTransmissibilities()
end

