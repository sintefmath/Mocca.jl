using Parameters

@with_kw struct PressurationBC
   trans::Float64 = 1.0
end

# function JutulDarcy.apply_flow_bc!(
#     acc,
#     q,
#     bc,
#     model::Jutul.SimulationModel{<:Any,T},
#     state,
#     time,
# ) where {T<:AdsorptionFlowSystem}
#     sys = model.system
#     mu = sys.p.fluid_viscosity
#     concentrations = state.concentrations
#     ctot = state.cTot
#     nph = length(acc)

#     rho_inj = bc.density
#     f_inj = bc.fractional_flow
#     c = bc.cell
#     mobility = 1 / mu

#     if q > 0
#         for component in eachindex(acc)
#             c_i = mobility * q * concentrations[component, c]
#             acc[component] += c_i
#         end
#     else
#         for component in eachindex(acc)
#             c_i = mobility * q * bc.fractional_flow[component]
#             acc[component] += c_i
#         end
#     end

# end

function pressure_left(state, model::AdsorptionFlowSystem, ::PressurationBC, t) 
    PH = model.p.p_high
    PL = model.p.p_low
    λ = model.p.λ
    # TODO: Shouldn't this be multiplied by 1/P0 
    return (PH - (PH - PL)*exp(-λ*t))
end

function pressure_right(state, model::AdsorptionFlowSystem, ::PressurationBC, t) 
    # Neumann
    return state.Pressure[end]
end

function mole_fraction_left(state, model::AdsorptionFlowSystem, ::PressurationBC) 
    return model.p.y_feed
end

function temperature_left(state, model::AdsorptionFlowSystem, ::PressurationBC) 
    return model.p.T_feed
end

function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::PressurationBC,
    time,
)

    # TODO: Refactor this a bit. Should probably have one function per side?

    @info "In adsorption bc" force
    state = storage.state
    #for bc in force
    sys = model.system
    
    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmisibility = force.trans

    # left side
    begin
        cell_left = 1
        P_left = pressure_left(state, sys, force, time)
        P = state.Pressure[cell_left]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P-P_left)
        y_left = mole_fraction_left(state, sys, force)
        y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)

        acc_i = view(acc, :, cell_left)
        cTot = P / (T_left * sys.p.R)
        c = y .* cTot
        acc_i[:] .+= cTot .* q .* (y_left .- y) .+ q .* c
    end
    # right side
    begin
        cell_right = 1
        P_right = pressure_right(state, sys, force, time)
        P = state.Pressure[cell_right]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P_right - P)
        y_left = mole_fraction_left(state, sys, force)
        y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)

        acc_i = view(acc, :, cell_right)
        cTot = P / (T_left * sys.p.R)
        c = y .* cTot
        acc_i[:] .+= cTot .* q .* (y_left .- y) .+ q .* c
    end
end

# function Jutul.apply_forces_to_equation!(
#     acc,
#     storage,
#     model::AdsorptionFlowModel,
#     eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
#     eq_s,
#     force::PressurationBC,
#     time,
# )
#     @info "In column bc"
#     # state = storage.state
#     # p = state.Pressure
#     # for bc in force
#     #     c = bc.cell
#     #     T_f = bc.trans_flow
#     #     Δp = p[c] - bc.pressure
#     #     q = T_f*Δp
#     #     acc_i = view(acc, :, c)
#     #     apply_flow_bc!(acc_i, q, bc, model, state, time)
#     # end
# end

# function Jutul.apply_forces_to_equation!(
#     acc,
#     storage,
#     model::AdsorptionFlowModel,
#     eq::Jutul.ConservationLaw{:WallConservedEnergy},
#     eq_s,
#     force::PressurationBC,
#     time,
# )
#     @info "In wall bc"
#     # state = storage.state
#     # p = state.Pressure
#     # for bc in force
#     #     c = bc.cell
#     #     T_f = bc.trans_flow
#     #     Δp = p[c] - bc.pressure
#     #     q = T_f*Δp
#     #     acc_i = view(acc, :, c)
#     #     apply_flow_bc!(acc_i, q, bc, model, state, time)
#     # end
# end