using Parameters

@with_kw struct PressurationBC
   trans::Float64 = 1.0
end

pressure_function(ph, pl, λ, t) = (ph - (ph-pl))*exp(-λ * t)

function pressure_left(state, model::AdsorptionFlowSystem, ::PressurationBC, t) 
    PH = model.p.p_high
    PL = model.p.p_low
    λ = model.p.λ
    # TODO: Shouldn't this be multiplied by 1/P0 
    return pressure_function(PH, PL, λ, t)
end

function pressure_right(state, model::AdsorptionFlowSystem, ::PressurationBC, t) 
    # TODO: This is how it's done in the matlab code at the moment, but this can't be right??
    PH = model.p.p_high
    PL = model.p.p_low
    λ = model.p.λ
    # TODO: Shouldn't this be multiplied by 1/P0 
    return pressure_function(PH, PL, λ, t)
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
    state = storage.state
    #for bc in force
    sys = model.system
    
    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmisibility = force.trans

    # left side
    begin
        cell_left = 1
        Δx = compute_dx(model, cell_left)

        P_left = pressure_left(state, sys, force, time)
        P = state.Pressure[cell_left]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P-P_left)
        y_left = mole_fraction_left(state, sys, force)
        y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)

        acc_i = view(acc, :, cell_left)
        cTot = P_left / (T_left * sys.p.R)
        c = y_left .* cTot
        acc_i[:] .+= (cTot .* q .* (y_left .- y) .+ q .* c)
    end
    # right side
    begin
        cell_right = size(state.Pressure, 1)
        Δx = compute_dx(model, cell_right)

        #@info "Flow" cell_right

        P_right = pressure_right(state, sys, force, time)
        P = state.Pressure[cell_right]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P_right - P)
        
        acc_i = view(acc, :, cell_right)
        c = state.concentrations[:, cell_right]
        #acc_i[:] .+= q .* c 
    end
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

    # TODO: Refactor this a bit. Should probably have one function per side? Can also be 
    # combined with the BC for the other equations.

    state = storage.state
    #for bc in force
    sys = model.system
    
    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmisibility = force.trans
    R = sys.p.R
    

    # left side
    begin
        ρ_g = sys.p.ρ_s
        cell_left = 1
        P_left = pressure_left(state, sys, force, time)

        Δx = compute_dx(model, cell_left)

        P = state.Pressure[cell_left]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P-P_left)
        #y = state.y[:, cell_left]
        T_left = temperature_left(state, sys, force)
        T = state.Temperature[cell_left]

        acc_i = view(acc, :, cell_left)
        #cTot = P / (T_left * sys.p.R)
        #c = y .* cTot
        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]
        
        acc_i[:] .+= ((q.*ρ_g.*C_pg.*(T_left - T)) + (q.*P_left./(R)).*C_pg.*avm)
    end
    # right side
    begin
        cell_right = size(state.Temperature, 1)
        Δx = compute_dx(model, cell_right)

      #  @info "Temp" cell_right

        P_right = pressure_right(state, sys, force, time)
        P = state.Pressure[cell_right]
        # v = -(T_{ij}/μ) ∇p
        q = -transmisibility * mobility * (P_right - P)
        
        acc_i = view(acc, :, cell_right)
        C_pg = state.C_pg[cell_right]
        avm = state.avm[cell_right]
        
        #acc_i[:] .+= (q.*P./(R).*C_pg.*avm)
    end
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

    # TODO: Refactor this a bit. Should probably have one function per side? Can also be 
    # combined with the BC for the other equations.

    state = storage.state
    #for bc in force
    sys = model.system
    
    μ = sys.p.fluid_viscosity
    mobility = 1.0 / μ
    transmisibility = force.trans
    R = sys.p.R

    # left side
    
    begin
        
        cell_left = 1
        Δx = compute_dx(model, cell_left)
        acc_i = view(acc, :, cell_left)
        T_left = temperature_left(state, sys, force)
        T = state.WallTemperature[cell_left]
        K_w = sys.p.K_w

        acc_i[:] .+= K_w * (T - T_left) / Δx
    end
    # right side
    begin
        cell_right = size(state.WallTemperature, 1)
        Δx = compute_dx(model, cell_right)
        #@info "Wall" cell_right
        acc_i = view(acc, :, cell_right)
        # FIXME: This might change in the future.
        T_right = temperature_left(state, sys, force) # Same temperature
        T = state.WallTemperature[cell_right]
        K_w = sys.p.K_w


        #acc_i[:] .-= K_w * (T_right - T) / Δx
    end
end
