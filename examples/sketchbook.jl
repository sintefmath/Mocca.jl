
module Mocca
import Jutul
import JutulDarcy

struct AdsorptionFlowSystem <: JutulDarcy.JutulSystem
end


JutulDarcy.number_of_components(sys::AdsorptionFlowSystem) = 2
JutulDarcy.has_other_phase(::AdsorptionFlowSystem) = false

struct BetterThermalSystem <: Jutul.JutulSystem
end

function Jutul.select_primary_variables!(S, ::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    S[:Pressure] = JutulDarcy.Pressure()
    S[:yCO2] = JutulDarcy.OverallMoleFractions()
    S[:qCO2] = JutulDarcy.TotalMass()
    S[:qN2] = JutulDarcy.TotalMass()
end

function Jutul.select_secondary_variables!(S, ::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    S[:cTot] = JutulDarcy.TotalMass()
    S[:cCO2] = JutulDarcy.TotalMass()
    S[:cN2] = JutulDarcy.TotalMass()
    S[:avm] = JutulDarcy.TotalMass()
    
end

function Jutul.select_equations!(eqs, sys::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    fdisc = model.domain.discretizations.mass_flow
    nc = JutulDarcy.number_of_components(sys)
    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :TotalMasses, nc)
end

function Jutul.select_parameters!(S, ::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    S[:Temperature] = JutulDarcy.Temperature()
    S[:fluidVolume] = JutulDarcy.FluidVolume()
    S[:bulkVolume] = JutulDarcy.BulkVolume() # solidVolume = buldVolume - fluidVolume

    # TODO: Find better type for Dispersion
    S[:axialDispersion] = JutulDarcy.Pressure()

    # TODO: Find proper type for fluidViscosity
    S[:fluidViscosity] = JutulDarcy.Transmissibilities()

end
time = 1.0
nc = 10
nstep=8
general_ad = false
T = time
tstep = repeat([T/nstep], nstep)
G = JutulDarcy.get_1d_reservoir(nc, general_ad = general_ad)
sys = AdsorptionFlowSystem()

model = Jutul.SimulationModel(G, sys)
@show model
end