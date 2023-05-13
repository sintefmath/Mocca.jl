struct PressurationBC


end

function JutulDarcy.apply_flow_bc!(
    acc,
    q,
    bc,
    model::Jutul.SimulationModel{<:Any,T},
    state,
    time,
) where {T<:AdsorptionFlowSystem}
    sys = model.system
    mu = sys.p.fluid_viscosity
    concentrations = state.concentrations
    ctot = state.cTot
    nph = length(acc)

    rho_inj = bc.density
    f_inj = bc.fractional_flow
    c = bc.cell
    mobility = 1 / mu

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


function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:adsorptionRates},
    eq_s,
    force::PressurationBC,
    time,
)
    @info "In adsorption bc"
    # state = storage.state
    # p = state.Pressure
    # for bc in force
    #     c = bc.cell
    #     T_f = bc.trans_flow
    #     Δp = p[c] - bc.pressure
    #     q = T_f*Δp
    #     acc_i = view(acc, :, c)
    #     apply_flow_bc!(acc_i, q, bc, model, state, time)
    # end
end

function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::PressurationBC,
    time,
)
    @info "In column bc"
    # state = storage.state
    # p = state.Pressure
    # for bc in force
    #     c = bc.cell
    #     T_f = bc.trans_flow
    #     Δp = p[c] - bc.pressure
    #     q = T_f*Δp
    #     acc_i = view(acc, :, c)
    #     apply_flow_bc!(acc_i, q, bc, model, state, time)
    # end
end

function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::PressurationBC,
    time,
)
    @info "In wall bc"
    # state = storage.state
    # p = state.Pressure
    # for bc in force
    #     c = bc.cell
    #     T_f = bc.trans_flow
    #     Δp = p[c] - bc.pressure
    #     q = T_f*Δp
    #     acc_i = view(acc, :, c)
    #     apply_flow_bc!(acc_i, q, bc, model, state, time)
    # end
end