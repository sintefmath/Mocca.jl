
module Mocca
import Jutul
import JutulDarcy
using StaticArrays
using Parameters
using Tullio

@with_kw struct AdsorptionFlowSystem <: JutulDarcy.MultiComponentSystem
    number_of_components::Int64 = 2
    molecularMassOfCO2::Float64 = 44.01e-3 # kg / mole
    molecularMassOfN2::Float64 = 28e-3 # kg/mole
    R::Float64 = 8.3144598 # J⋅mol^−1⋅K^−1.
    rho_ref::Float64 = 1.0 # TODO: DELETEME!!!
    Φ::Float64 = 0.37 # TODO: We should not hardcode this....
end

struct AdsorptionRates <: Jutul.VectorVariables
end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AdsorptionRates) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AdsorptionRates) = JutulDarcy.number_of_components(model.system)

struct AverageMolecularMass <: Jutul.ScalarVariable
end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AverageMolecularMass) = 1

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AverageMolecularMass) = 1

struct Concentrations <: Jutul.VectorVariables end


Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

JutulDarcy.phase_names(::AdsorptionFlowSystem) = ["ouronlyphase"]
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
    S[:Pressure] = JutulDarcy.Pressure(minimum=π) # FIXME: Proper lower value 
    S[:y] = JutulDarcy.OverallMoleFractions()
end

function Jutul.select_secondary_variables!(
    S,
    system::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    S[:adsorptionRates] = AdsorptionRates()
    #S[:qCO2] = JutulDarcy.TotalMass()
    #S[:qN2] = JutulDarcy.TotalMass()
    S[:cTot] = JutulDarcy.TotalMass()
    S[:concentrations] = Concentrations()
    S[:avm] = AverageMolecularMass()
    S[:TotalMasses] = JutulDarcy.TotalMasses()

    # 🙏
    # Hope this works... Should provide uniqueness for system.
    nph = JutulDarcy.number_of_phases(system)
    S[:PhaseMassDensities] = JutulDarcy.ConstantCompressibilityDensities(nph)

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
    S[:molecularMassOfCO2] = JutulDarcy.TotalMass()
    S[:molecularMassOfN2] = JutulDarcy.TotalMass()
    S[:PhaseViscosities] = JutulDarcy.PhaseViscosities()

end

const CO2INDEX = 1
const N2INDEX = 2

const AdsorptionFlowModel = Jutul.SimulationModel{<:Any, <:AdsorptionFlowSystem, <:Any, <:Any}

# We need this since the jutul_secondary macro assumes
# update_secondary_variable! is in the namespace
# import Jutul: update_secondary_variable!


Jutul.@jutul_secondary function update_our_total_masses!(
    totmass,
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    y,
    cTot,
    PhaseMassDensities, 
    ix
)
    println("Updating total mass")
    sys = model.system

    for cell in ix
        #@info "Data in cell" cTot[cell] y[:, cell] PhaseMassDensities[1, cell]  model.domain.grid.pore_volumes[cell]
        totmass[1, cell] = cTot[cell] * y[1, cell] * PhaseMassDensities[1, cell] * model.domain.grid.pore_volumes[cell]
        totmass[2, cell] = cTot[cell] * y[2, cell] * PhaseMassDensities[1, cell] * model.domain.grid.pore_volumes[cell]
        # totmass[CO2COMPONENTINDEX, cell] = y
        # totmass[N2COMPONENTINDEX, cell] = y
    end
end

function JutulDarcy.component_mass_fluxes!(q, face, state, model::Jutul.SimulationModel{G, S}, kgrad, upw) where {G<:Any, S<:AdsorptionFlowSystem}
    # This is defined for us:
    # kgrad = TPFA(left, right, face_sign)
    # upw = SPU(left, right)
    # TODO: Implement me


    disc = JutulDarcy.kgrad_common(face, state, model, kgrad)
    (∇p, T_f, gΔz) = disc
    # @show T_f
    # @show ∇p
    for ph in eachindex(q)
        q_darcy = -T_f * ∇p
        
        q = setindex(q, q_darcy, ph)
    end
    # @info "Flux $face" q
    return q
end

# @inline function Jutul.face_flux!(Q, left, right, face, face_sign, eq::JutulDarcy.ConservationLaw{:TotalMasses}, state, model::AdsorptionFlowModel, dt, flow_disc::Jutul.TwoPointPotentialFlowHardCoded)
#     # # Specific version for tpfa flux
#     # # TODO: Add general version for thermal
#     # grad = TPFA(left, right, face_sign)
#     # upw = SPU(left, right)
#     # q = thermal_heat_flux!(face, state, model, grad, upw)
#     # return setindex(Q, q, 1)
#     println("Hello")
#     return Q
# end

function Jutul.convergence_criterion(model::Jutul.SimulationModel{D, S}, storage, eq::Jutul.ConservationLaw, eq_s, r; dt = 1) where {D, S<:AdsorptionFlowSystem}
    M = Jutul.global_map(model.domain)
    
    v = x -> Jutul.as_value(Jutul.active_view(x, M, for_variables = false))
    #Φ = v(storage.state.FluidVolume)
    Φ = model.system.Φ
    
    ρ = v(storage.state.PhaseMassDensities)

    @show size(r)
    @show size(ρ)
    e = [maximum(abs.(r[j, :]) * dt / (ρ[1, :]*Φ)) for j in 1:2]
    
    N = length(Φ)
    pv_t = sum(Φ)
    avg_density = sum(ρ, dims = 2)./N
    r_sum = sum(r, dims = 2)
    mb = @. (dt/pv_t)*abs(r_sum)/avg_density

    names = ["CO2", "N2"]
    R = (CNV = (errors = e, names = names),
         MB = (errors = mb, names = names))
    return R
end



Jutul.@jutul_secondary function update_adsorption_rates!(totmass, tv::AdsorptionRates,
     model::Jutul.SimulationModel{G, S}, ix) where {G, S<:AdsorptionFlowSystem}
    # TODO: Implement this
    println("Updating adsorption rates")
end

Jutul.@jutul_secondary function update_cTot!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G, S}, Pressure, Temperature,  ix) where {G, S<:AdsorptionFlowSystem}
    # Update cTot
    sys = model.system

    for cellindex in ix
        ctot[cellindex] = Pressure[cellindex] / (sys.R * Temperature[cellindex])
        #@info "Updating cTot to " ctot[cellindex] Pressure[cellindex] sys.R Temperature[cellindex]

    end
end

Jutul.@jutul_secondary function update_avm!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S<:AdsorptionFlowSystem}
    println("Updating avm")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        sys = model.system
        molecularMassOfCO2 = sys.molecularMassOfCO2
        molecularMassOfN2 = sys.molecularMassOfN2
        avm[cellindex] = y[CO2INDEX, cellindex] * molecularMassOfCO2 + y[N2INDEX, cellindex] * molecularMassOfN2
    end
end

Jutul.@jutul_secondary function update_concentrations!(concentrations, tv::Concentrations, model::Jutul.SimulationModel{G, S}, y, cTot, ix) where {G, S<:AdsorptionFlowSystem}
    println("Updating concentrations")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        for component in 1:JutulDarcy.number_of_components(model.system)
            concentrations[component, cellindex] = y[component, cellindex] * cTot[cellindex]
        end
    end
end

# include("dontdothis.jl")




time = 1.0
nc = 10
nstep = 8
general_ad = false
T = time
tstep = repeat([T / nstep], nstep)
timesteps = tstep*3600*24 # Convert time-steps from days to seconds

ϵ = 0.37
r_in = 0.289/2.0
perm = 4/150 * (ϵ/(1-ϵ))^2*r_in^2

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
irate = 500*sum(G.grid.pore_volumes)/time
# src  = [JutulDarcy.SourceTerm(1, irate, fractional_flow = [1.0, 0.0]), 
#     JutulDarcy.SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
# forces = JutulDarcy.setup_forces(model, sources = src)
@info "parameter set" parameters model.domain.grid.trans
state0 = Jutul.setup_state(model, Pressure = p0, y = [0.0, 1.0])
# Simulate and return
sim = Jutul.Simulator(model, state0 = state0, parameters = parameters)
states, report = Jutul.simulate(sim, timesteps, info_level = 5)#, forces=forces)
end

Mocca.model