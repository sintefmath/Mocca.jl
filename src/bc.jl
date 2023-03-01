
function JutulDarcy.apply_flow_bc!(acc, q, bc, model::Jutul.SimulationModel{<:Any,T}, state, time) where {T<:AdsorptionFlowSystem}

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
    mobility = 1 / mu[1, c]

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