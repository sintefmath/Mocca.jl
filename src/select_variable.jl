function Jutul.select_primary_variables!(
    S,
    system::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    S[:Pressure] = JutulDarcy.Pressure(minimum=π) # FIXME: Proper lower value 
    S[:y] = GasMoleFractions()
    S[:adsorptionRates] = AdsorptionRates()
    S[:Temperature] = JutulDarcy.Temperature()
    S[:WallTemperature] = JutulDarcy.Temperature()
end

function Jutul.select_secondary_variables!(
    S,
    system::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    S[:cTot] = JutulDarcy.TotalMass()
    S[:concentrations] = Concentrations()
    # Not using AVM at the moment, so disabling. 
    # TODO: Renable AVM
    S[:avm] = AverageMolecularMass()
    S[:TotalMasses] = JutulDarcy.TotalMasses()
    S[:AdsorptionMassTransfer] = AdsorptionMassTransfer()
    # 🙏 # Might still need this.
    # Hope this works... Should provide uniqueness for system. We don't seem to need this anymore
    #nph = JutulDarcy.number_of_phases(system)
    #S[:PhaseMassDensities] = JutulDarcy.ConstantCompressibilityDensities(nph)

    # For the energy equations
    S[:ColumnConservedEnergy] = ColumnEnergy()
    S[:WallConservedEnergy] = WallEnergy()
    S[:ΔH] = EnthalpyChange()
    S[:C_pa] = SpecificHeatCapasityAdsorbent()
    S[:C_pg] = SpecificHeatCapasityFluid()
end

function Jutul.select_equations!(
    eqs,
    sys::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    fdisc = model.domain.discretizations.mass_flow
    nc = JutulDarcy.number_of_components(sys)
    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :TotalMasses, nc)
    eqs[:mass_transfer] = Jutul.ConservationLaw(fdisc, :adsorptionRates, nc)

    eqs[:energy_column] = Jutul.ConservationLaw(fdisc, :ColumnConservedEnergy, 1)
    eqs[:energy_wall] = Jutul.ConservationLaw(fdisc, :WallConservedEnergy, 1)
end

function Jutul.default_value(model::AdsorptionFlowModel, ::JutulDarcy.BulkVolume)
    Φ = model.system.p.Φ
    error()
end
function Jutul.select_parameters!(S, ::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    S[:solidVolume] = JutulDarcy.BulkVolume()
    S[:fluidVolume] = JutulDarcy.FluidVolume()

end

