
module Mocca
import Jutul
import JutulDarcy

using Parameters

@with_kw struct AdsorptionFlowSystem <: JutulDarcy.MultiComponentSystem
    number_of_components::Int64 = 2
end

JutulDarcy.number_of_phases(::AdsorptionFlowSystem) = 1

# function JutulDarcy.degrees_of_freedom_per_entity(
#     model::Jutul.SimulationModel{G,S},
#     v::JutulDarcy.TotalMasses,
# ) where {G<:Any,S<:AdsorptionFlowSystem}
#     JutulDarcy.number_of_phases(model.system)
# end


JutulDarcy.number_of_components(sys::AdsorptionFlowSystem) = sys.number_of_components
JutulDarcy.has_other_phase(::AdsorptionFlowSystem) = false

function Jutul.select_primary_variables!(
    S,
    ::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    S[:Pressure] = JutulDarcy.Pressure()
    S[:y] = JutulDarcy.OverallMoleFractions()
end

function Jutul.select_secondary_variables!(
    S,
    ::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    S[:adsorptionRateCO2] = JutulDarcy.TotalMass() # TODO: Make MassPerSecondType
    S[:adsorptionRateN2] = JutulDarcy.TotalMass() # TODO: Make MassPerSecondType
    #S[:qCO2] = JutulDarcy.TotalMass()
    #S[:qN2] = JutulDarcy.TotalMass()
    S[:cTot] = JutulDarcy.TotalMass()
    S[:cCO2] = JutulDarcy.TotalMass()
    S[:cN2] = JutulDarcy.TotalMass()
    S[:avm] = JutulDarcy.TotalMass()
    S[:TotalMasses] = JutulDarcy.TotalMasses()
end

function Jutul.select_equations!(
    eqs,
    sys::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    fdisc = model.domain.discretizations.mass_flow
    nc = JutulDarcy.number_of_components(sys)
    eqs[:mass_conservation] = Jutul.ConservationLaw(fdisc, :TotalMasses, nc)
end

function Jutul.select_parameters!(S, ::AdsorptionFlowSystem, model::Jutul.SimulationModel)
    S[:Temperature] = JutulDarcy.Temperature()

    # TODO: Find better type for Dispersion
    S[:axialDispersion] = JutulDarcy.Pressure()

    # TODO: Find proper type for fluidViscosity
    S[:fluidViscosity] = JutulDarcy.Transmissibilities()
    S[:solidVolume] = JutulDarcy.BulkVolume()

end

const CO2COMPONENTINDEX = 1
const N2COMPONENTINDEX = 2

import Jutul: update_secondary_variable!
Jutul.@jutul_secondary function update_our_total_masses!(
    totmass,
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    adsorptionRateCO2,
    adsorptionRateN2,
    y,
    cTot,
    cCO2,
    cN2,
    axialDispersion,
    ix
)
    sys = model.system
    # for cell in ix
    #     totmass[CO2COMPONENTINDEX, cell] = y
    #     totmass[N2COMPONENTINDEX, cell] = y
    # end
end

Jutul.@jutul_secondary function update_total_masses!(totmass, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G, S}, ix) where {G, S<:AdsorptionFlowSystem}
    error("Implement me!!!")
end

time = 1.0
nc = 10
nstep = 8
general_ad = false
T = time
tstep = repeat([T / nstep], nstep)
timesteps = tstep*3600*24 # Convert time-steps from days to seconds

ϵ = 0.37
r_in = 0.289/2.0
perm = -4/150 * (ϵ/(1-ϵ))^2*r_in^2

G = JutulDarcy.get_1d_reservoir(nc, general_ad = general_ad, poro=ϵ , perm=perm)
sys = AdsorptionFlowSystem()

model = Jutul.SimulationModel(G, sys)

# pv = vol * poro
pv = model.domain.grid.pore_volumes
volumes = pv / ϵ 
solid_volume = volumes * (1 - ϵ )

# axial Dispersion
dp = 0.002
Dm = 1.6e-5
V0_inter = 0.03653;         # Interstitial inlet velocity [m/s]
V0 = V0_inter*ϵ;         # Inlet velocity [m/s]

DL = 0.7 * Dm + 0.5 * V0 * dp
bar = 1e5
p0 = 100*bar
parameters = Jutul.setup_parameters(model, Temperature=298,
    solidVolume = solid_volume, 
    axialDispersion = DL,
    fluidViscosity = 1.72e-5)
state0 = Jutul.setup_state(model, Pressure = p0, y = [0.0, 1.0])
# Simulate and return
sim = Jutul.Simulator(model, state0 = state0, parameters = parameters)
#states, report = Jutul.simulate(sim, timesteps)
end

Mocca.model