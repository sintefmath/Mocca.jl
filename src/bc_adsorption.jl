@with_kw struct AdsorptionBC
    trans::Float64 = 1.0
end



function velocity_left(state, model::AdsorptionFlowSystem, ::PressurationBC, t)
    return model.p.v_feed
end

function mole_fraction_left(state, model::AdsorptionFlowSystem, ::PressurationBC)
    return model.p.y_feed
end

function temperature_left(state, model::AdsorptionFlowSystem, ::PressurationBC)
    return model.p.T_feed
end

function pressure_right(state, model::AdsorptionFlowSystem, ::PressurationBC, t)
    return  model.p.p_high
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::AdsorptionBC,
    time,
)

    # TODO: Refactor this a bit. Should probably have one function per side?
    state = storage.state
    #for bc in force
    sys = model.system

    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmissibility = force.trans

    # left side
    begin
        cell_left = 1
        Δx = compute_dx(model, cell_left)

        q = - velocity_left * A_f
        P = state.Pressure[cell_left]
        P_left = q / tranmissibility / mobility + P

        y_left = mole_fraction_left(state, sys, force)
        y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)

        acc_i = view(acc, :, cell_left)
        cTot = P_left / (T_left * sys.p.R)
        c = y_left .* cTot
        for i in eachindex(y)

            mysource = cTot*q*(y_left[i] - y[i]) + q*c[i]
            acc[i, cell_left] -= mysource
        end
    end

    # left side
    begin
        cell_right = 30 # TODO: Don't hardcode final index!
        Δx = compute_dx(model, cell_right)

        P_right = pressure_right(state, sys, force, time)
        P = state.Pressure[cell_right]
        T = state.Temperature[cell_right]
        cTot = P / (T * sys.p.R)

        q = -transmissibility * mobility * (P_right - P)

        for i in eachindex(y)

            mysource =  q*c[i]
            acc[i, cell_left] -= mysource
        end
    end
  
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::AdsorptionBC,
    time,
)

    # TODO: Refactor this a bit. Should probably have one function per side? Can also be 
    # combined with the BC for the other equations.

    state = storage.state
    #for bc in force
    sys = model.system

    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmissibility = force.trans
    R = sys.p.R


    # left side
    begin
        ρ_g = sys.p.ρ_s
        cell_left = 1
        P_left = pressure_left(state, sys, force, time)

        P = state.Pressure[cell_left]
        # v = -(T_{ij}/μ) ∇p
        q = -transmissibility * mobility * (P - P_left)
        #y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)
        T = state.Temperature[cell_left]

        acc_i = view(acc, :, cell_left)
        #cTot = P / (T_left * sys.p.R)
        #c = y .* cTot
        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]


        # @info "BC " q (T_left - T) ((q .* ρ_g .* C_pg .* (T_left - T)) + (q .* P_left ./ (R)) .* C_pg .* avm) acc_i[:]
        acc[cell_left] += ((q * ρ_g * C_pg * (T_left - T)) + (q * P_left / (R)) * C_pg * avm)
    end
   
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::AdsorptionBC,
    time,
)

    # TODO: Refactor this a bit. Should probably have one function per side? Can also be 
    # combined with the BC for the other equations.

    state = storage.state
    #for bc in force
    sys = model.system

    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmissibility = force.trans
    R = sys.p.R

    # left side

    begin

        cell_left = 1
        Δx = compute_dx(model, cell_left)
        # acc_i = view(acc, :, cell_left)
        T_left = sys.p.T_a
        T = state.WallTemperature[cell_left]
        K_w = sys.p.K_w

        A_w = area_wall(sys)
        acc[cell_left] -= A_w * K_w * (T - T_left) / Δx
    end
  
end