using Parameters

@with_kw struct PressurisationBC{T, N}
    y_feed::SVector{N,T}
    PH::T
    PL::T
    λ::T
    T_feed::T
    cell_left::Int
end

function pressure_left(force::PressurisationBC, time)
    PH = force.PH
    PL = force.PL
    λ = force.λ
    return (PH - (PH - PL) * exp(-λ * time))
end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::PressurisationBC,
    time,
)

    state = storage.state

    pars = model.system.p
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model)

    # left side
    begin
        cell_left = force.cell_left
        P = state.Pressure[cell_left]
        y = state.y[:, cell_left] 

        P_bc = pressure_left(force, time)
        y_bc = force.y_feed        
        T_bc = force.T_feed

        q = -trans * mob * (P - P_bc)

        cTot = P_bc / (T_bc * R)
        c = y_bc .* cTot

        for i in eachindex(y)
            mysource = cTot * q * (y_bc[i] - y[i]) + q * c[i]
            acc[i, cell_left] += mysource
        end

    end
  
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::PressurisationBC,
    time,
)

    state = storage.state


    pars = model.system.p
    ρ_g = pars.ρ_g
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model)

    # left side
    begin
        cell_left = force.cell_left
        P = state.Pressure[cell_left]
        y = state.y[:, cell_left] 
        T = state.Temperature[cell_left]
        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]        

        P_bc = pressure_left(force, time)
        y_bc = force.y_feed        
        T_bc = force.T_feed

        q = -trans * mob * (P - P_bc)

        cTot = P_bc / (T_bc * pars.R)
        c = y_bc .* cTot

        bc_src = -((q * ρ_g * C_pg * (T_bc - T)) + (q * P_bc / R) * C_pg * avm)
        acc[cell_left] -= bc_src

    end

end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::PressurisationBC,
    time,
)

    state = storage.state

    pars = model.system.p
 
    # left side
    begin
        cell_left = force.cell_left
        trans_wall = calc_bc_wall_trans(model)

        T = state.WallTemperature[cell_left]
        T_bc = pars.T_a

        bc_src = -(trans_wall * (T - T_bc))
        acc[cell_left] -= bc_src
    end
end
