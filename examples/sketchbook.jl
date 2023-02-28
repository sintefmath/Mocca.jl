
module Mocca
import Jutul
import JutulDarcy
using StaticArrays
using Parameters
using Tullio
using CairoMakie
using MakiePublication
using LinearAlgebra
@with_kw struct AdsorptionFlowSystem <: JutulDarcy.MultiComponentSystem
    number_of_components::Int64 = 2
    molecularMassOfCO2::Float64 = 44.01e-3 # kg / mole
    molecularMassOfN2::Float64 = 28e-3 # kg/mole
    R::Float64 = 8.3144598 # J⋅mol^−1⋅K^−1.
    rho_ref::Float64 = 1.0 # TODO: DELETEME!!!
    Φ::Float64 = 0.37 # TODO: We should not hardcode this....
    b0::SVector{2, Float64} = @SVector [8.65e-7, 2.5e-6]
    d0::SVector{2, Float64} = @SVector [2.63e-8, 0.0]
    ΔUbi::SVector{2, Float64} = @SVector [-36_641.21, -1.58e4]
    ΔUdi::SVector{2, Float64} = @SVector [-35_690.66, 0.0]
    qsbi::SVector{2, Float64} = @SVector [3489.44, 6613.551]
    qsdi::SVector{2, Float64} = @SVector [2872.35, 0.00]
    ϵ_p::Float64 = 0.35
    D_m::Float64 = 1.6e-5
    τ::Float64 = 3.0
    d_p::Float64 = 0.002
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

struct AdsorptionMassTransfer <: Jutul.VectorVariables end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

JutulDarcy.phase_names(::AdsorptionFlowSystem) = ["ouronlyphase"]
JutulDarcy.number_of_phases(::AdsorptionFlowSystem) = 1

struct GasMoleFractions <: JutulDarcy.CompositionalFractions
    dz_max::Float64
    GasMoleFractions(;dz_max = 0.2) = new(dz_max)
end

function Jutul.minimum_value(::GasMoleFractions)
    return 1e-20
end

# Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::GasMoleFractions) = JutulDarcy.number_of_components(model.system)

# Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::GasMoleFractions) = JutulDarcy.number_of_components(model.system)

Jutul.absolute_increment_limit(z::GasMoleFractions) = z.dz_max


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
    S[:y] = GasMoleFractions()
    S[:adsorptionRates] = AdsorptionRates()

end

function Jutul.select_secondary_variables!(
    S,
    system::AdsorptionFlowSystem,
    model::Jutul.SimulationModel,
)
    #S[:qCO2] = JutulDarcy.TotalMass()
    #S[:qN2] = JutulDarcy.TotalMass()
    S[:cTot] = JutulDarcy.TotalMass()
    S[:concentrations] = Concentrations()
    # Not using AVM at the moment, so disabling. 
    # TODO: Renable AVM
    #S[:avm] = AverageMolecularMass()
    S[:TotalMasses] = JutulDarcy.TotalMasses()
    S[:AdsorptionMassTransfer] = AdsorptionMassTransfer()
    # 🙏 # Might still need this.
    # Hope this works... Should provide uniqueness for system. We don't seem to need this anymore
    #nph = JutulDarcy.number_of_phases(system)
    #S[:PhaseMassDensities] = JutulDarcy.ConstantCompressibilityDensities(nph)

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
    # S[:Pressure] = JutulDarcy.Pressure(minimum=π) # FIXME: Proper lower value 
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
    concentrations, 
    ix
)
    sys = model.system
    pv = Jutul.physical_representation(model.domain).pore_volumes
    pv = sys.Φ * ones(size(pv)) # FIXME: Remove this
    for cell in ix
        #@info "Data in cell" cTot[cell] y[:, cell] PhaseMassDensities[1, cell]  model.domain.grid.pore_volumes[cell]
        totmass[1, cell] = concentrations[1, cell] * pv[cell]
        totmass[2, cell] = concentrations[2, cell] * pv[cell]
        # totmass[CO2COMPONENTINDEX, cell] = y
        # totmass[N2COMPONENTINDEX, cell] = y
    end
end



Jutul.@jutul_secondary function update_adsorption_mass_transfer(
    adsorption_mass_transfer,
    tv::AdsorptionMassTransfer,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    concentrations,
    Temperature, 
    adsorptionRates,
    ix
)
    for cell in ix
        qstar = compute_equilibrium(model.system, concentrations[:, cell], Temperature[cell])
        k = compute_ki(model.system, concentrations[:, cell], qstar)
        force = k.*(qstar .- adsorptionRates[:, cell])

        # @info "cell $ix" qstar k concentrations adsorptionRates force

        adsorption_mass_transfer[:, cell] = force
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
    c = state.concentrations
    μ = state.PhaseViscosities
    q_darcy = -T_f * ∇p

    R = sys.R
    L = kgrad.left
    R = kgrad.right
    
    favg(X) = (X[L] + X[R])/2
    P = favg(state.Pressure)
    T = favg(state.Temperature)
    C = P/(R*T)

    D_l = favg(state.axialDispersion)
    for component in eachindex(q)
        F_c = cell -> c[component, cell]/μ[cell]
        c_face = JutulDarcy.upwind(upw, F_c, q_darcy)
        y_i = view(state.y, component, :)
        q_i = c_face*q_darcy + C*D_l*JutulDarcy.gradient(y_i, kgrad)
        
        q = setindex(q, q_i, component)
    end
    # @info "Flux $face" q
    return q
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::Jutul.ConservationLaw{:TotalMasses}, model::AdsorptionFlowModel, Δt, ldisc = Jutul.local_discretization(eq, self_cell)) where T_e
    #@info "Updating total masses"
    
    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    # Compute ∇⋅V
    disc = eq.flow_discretization
    flux(face) = Jutul.face_flux(face, eq, state, model, Δt, disc, ldisc, Val(T_e))
    div_v = ldisc.div(flux)
    # @info "Forcing term" forcing_term
    for i in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, i, self_cell)
        # @info i ∂M∂t forcing_term[i, self_cell]
        eq_buf[i] = ∂M∂t + div_v[i] + forcing_term[i, self_cell]
    end
end

function Jutul.update_equation_in_entity!(eq_buf::AbstractVector{T_e}, self_cell, state, state0, eq::Jutul.ConservationLaw{:adsorptionRates}, model::AdsorptionFlowModel, Δt, ldisc = Jutul.local_discretization(eq, self_cell)) where T_e
    
    # Compute accumulation term
    conserved = Jutul.conserved_symbol(eq)
    M₀ = state0[conserved]
    M = state[conserved]

    forcing_term = state[:AdsorptionMassTransfer]
    ϵ = model.system.Φ
    # @info "From transfer system" forcing_term
    for component in eachindex(eq_buf)
        ∂M∂t = Jutul.accumulation_term(M, M₀, Δt, component, self_cell)
        eq_buf[component] = ∂M∂t - forcing_term[component, self_cell]/(1-ϵ)
    end
end

function compute_equilibrium(sys::AdsorptionFlowSystem, concentration, temperature)
    qstar = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    b = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    d = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    for i in 1:JutulDarcy.number_of_components(sys)
        b[i] = sys.b0[i] * exp(-sys.ΔUbi[i]/(sys.R*temperature))
        d[i] = sys.d0[i] * exp(-sys.ΔUdi[i]/(sys.R*temperature))
    end
    
    for i in 1:JutulDarcy.number_of_components(sys)
        qstar[i] = sys.qsbi[i]*b[i]*concentration[i] / (1 + sum(b.*concentration)) + sys.qsdi[i]*d[i]*concentration[i] / (1 + sum(d.*concentration))
    end
    return qstar
end

function compute_ki(sys::AdsorptionFlowSystem, concentration, qstar)
    D_p = sys.D_m / sys.τ
    r_p = sys.d_p / 2.0

    return concentration ./ qstar .* 15 .* sys.ϵ_p .* D_p ./ r_p^2
end

function Jutul.convergence_criterion(model::Jutul.SimulationModel{D, S}, storage, eq::Jutul.ConservationLaw, eq_s, r; dt = 1) where {D, S<:AdsorptionFlowSystem}
    M = Jutul.global_map(model.domain)
    
    v = x -> Jutul.as_value(Jutul.active_view(x, M, for_variables = false))
    #Φ = v(storage.state.FluidVolume)
    Φ = model.system.Φ
    
    scale = 1.0
    # scale = 1/1000.0
    # scale = 1/1e13
    e = scale.*[maximum(abs.(r[j, :]) * dt / (Φ)) for j in 1:2]

    # e = sum(abs, r, dims = 2)
    y = storage.state.y
    q = storage.state.adsorptionRates
    # @show r
    # @info "conv " size(r)
    #r_scaled = [r[:, cell] ./ [y[:, cell]..., q[:, cell]...] for cell in eachindex(q[1, :])]

    e  = norm(r)
    
    N = length(Φ)
    pv_t = sum(Φ)
    avg_density = 1.0#sum(ρ, dims = 2)./N
    r_sum = sum(r, dims = 2)
    # mb = @. (dt/pv_t)*abs(r_sum)/avg_density

    names = ["CO2", "N2"]
    R = (CNV = (errors = e, names = names),
         # MB = (errors = mb, names = names)
         )
    return R
end

Jutul.@jutul_secondary function update_cTot!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G, S}, Pressure, Temperature,  ix) where {G, S<:AdsorptionFlowSystem}
    # Update cTot
    sys = model.system

    for cellindex in ix
        ctot[cellindex] = Pressure[cellindex] / (sys.R * Temperature[cellindex])
    end
end

Jutul.@jutul_secondary function update_avm!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S<:AdsorptionFlowSystem}
    #println("Updating avm")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        sys = model.system
        molecularMassOfCO2 = sys.molecularMassOfCO2
        molecularMassOfN2 = sys.molecularMassOfN2
        avm[cellindex] = y[CO2INDEX, cellindex] * molecularMassOfCO2 + y[N2INDEX, cellindex] * molecularMassOfN2
    end
end

Jutul.@jutul_secondary function update_concentrations!(concentrations, tv::Concentrations, model::Jutul.SimulationModel{G, S}, y, cTot, ix) where {G, S<:AdsorptionFlowSystem}
   # println("Updating concentrations")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        for component in 1:JutulDarcy.number_of_components(model.system)
            concentrations[component, cellindex] = y[component, cellindex] * cTot[cellindex]
        end
    end
end

# include("dontdothis.jl")


function JutulDarcy.apply_flow_bc!(acc, q, bc, model::Jutul.SimulationModel{<:Any, T}, state, time) where T<:AdsorptionFlowSystem

    mu = state.PhaseViscosities
    # @show size(mu)
    # rho = state.PhaseMassDensities
    concentrations = state.concentrations
    ctot = state.cTot
    nph = length(acc)

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    
    # = 1/mu[1, c]
    mobility  = 1/mu[1, c]
    
    if q > 0
        for component in eachindex(acc)
            c_i = mobility * q * concentrations[component, c] 
            acc[component] += c_i
        end
    else
        for component in eachindex(acc)
            c_i = mobility * q * bc.fractional_flow[component]
            acc[component] += c_i
        end
    end
   
end

time = 1.0
nc = 10
nc = 2
nstep = 8
# nstep = 1000
general_ad = true
T = time
tstep = repeat([T / nstep], nstep)
# timesteps = tstep*3600*24 # Convert time-steps from days to seconds
timesteps = tstep

ϵ = 0.37
r_in = 0.289/2.0
perm = 4/150 * ((ϵ/(1-ϵ))^2)*r_in^2
perm = 1e-14

G = JutulDarcy.get_1d_reservoir(nc, general_ad = general_ad, poro=ϵ , perm=perm)
sys = AdsorptionFlowSystem()

ctx = Jutul.DefaultContext()
model = Jutul.SimulationModel(G, sys, context = ctx)

g = Jutul.physical_representation(model.domain)
# pv = vol * poro
pv = g.pore_volumes
volumes = pv / ϵ 
solid_volume = volumes * (1 - ϵ )

# axial Dispersion
dp = sys.d_p
Dm = sys.D_m
V0_inter = 0.03653;         # Interstitial inlet velocity [m/s]
V0 = V0_inter*ϵ;         # Inlet velocity [m/s]

DL = 0.7 * Dm + 0.5 * V0 * dp
bar = 1e5
p0 = 100*bar
parameters = Jutul.setup_parameters(model, Temperature=298,
    solidVolume = solid_volume, 
    axialDispersion = DL,
    fluidViscosity = 1.72e-5,
    PhaseMassDensities=1.0)

# TODO: Find a nicer way to specify trans on the boundary
d = JutulDarcy.FlowBoundaryCondition(nc, 2*p0, trans_flow = g.trans[1], fractional_flow = (0.5, 0.5))
# forces = Jutul.setup_forces(model, sources = [], bc = d)
forces = Jutul.setup_forces(model)
irate = 500*sum(g.pore_volumes)/time
# src  = [JutulDarcy.SourceTerm(1, irate, fractional_flow = [1.0, 0.0]), 
#     JutulDarcy.SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
# forces = JutulDarcy.setup_forces(model, sources = src)
@info "parameter set" parameters g.trans

yCO2 = 1e-15
# yCO2 = 1e-10
initY = [yCO2, 1 - yCO2]
ctot = p0/sys.R/298
c = ctot .* initY
equilinit = compute_equilibrium(sys, c, 298) # TODO: Should this still be zero for CO2?
equilinit = [0, equilinit[2]]
# equilinit = [0, 0.0]
@show equilinit


p_init = p0
p_init = [100*bar, 99.9999*bar]
state0 = Jutul.setup_state(model,
    Pressure = p_init, 
    y = initY, 
    adsorptionRates=equilinit)
# Simulate and return
sim = Jutul.Simulator(model, state0 = state0, parameters = parameters)
states, report = Jutul.simulate(sim, timesteps, info_level = 3, forces=forces, max_timestep_cuts = 0)
# states, report = Jutul.simulate(sim, timesteps, info_level = 4, forces=forces, max_timestep_cuts = 0, max_nonlinear_iterations = 0)

# with_theme(theme_web()) do 
#     f = CairoMakie.Figure()
#     ax = CairoMakie.Axis(f[1, 1], ylabel = "y", title = "Adsorption")
#     x = range(0, stop = 1, length = nc)
#     for i in 1:6:length(states)
#         CairoMakie.lines!(ax, x, states[i][:y][2, :], color = :darkgray)
#     end
#     display(f)
# end

with_theme(theme_web()) do 
    f = CairoMakie.Figure()
    ax = CairoMakie.Axis(f[1, 1], ylabel = "y", title = "Adsorption")
    y = map(x -> x[:y][1], states)
    for i in 1:6:length(states)
        CairoMakie.lines!(ax, y, color = :darkgray)
    end
    display(f)
end
display(Mocca.sim.storage.LinearizedSystem.jac)
display(Mocca.sim.storage.LinearizedSystem.jac[1:2:end,:])
end

Mocca.model