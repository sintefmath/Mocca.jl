using Parameters

@with_kw struct EvacuationBC{T}
    PL::T
    PI::T
    λ::T
    cell_left::Int
    cell_right::Int
    t_stage
end

function pressure_left(force::EvacuationBC, time)

    step_end = cumsum(force.t_stage)
    cycle_time = sum(force.t_stage)    
    cycle_no = floor(time/cycle_time)
 
    t_0 = cycle_no*cycle_time + step_end[3]
    t = time - t_0

    PL = force.PL
    PI = force.PI
    λ = force.λ

    return (PL + (PI - PL) * exp(-λ * t))
end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::EvacuationBC,
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
        T = state.Temperature[cell_left]
        y = state.y[:, cell_left] 

        P_bc = pressure_left(force, time)

        q = -trans * mob * (P_bc - P)

        cTot = P / (T * R)

        for i in eachindex(y)
            c = y[i] * cTot
            mysource =  -(q * c)
            acc[i, cell_left] -= mysource
        end

    end
  
end


function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::EvacuationBC,
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

        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]        

        P_bc = pressure_left(force, time)


        q = -trans * mob * (P_bc - P)


        bc_src = -(q * P / R * C_pg * avm)
        acc[cell_left] -= bc_src

    end

end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::EvacuationBC,
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

    # right side
    begin
        cell_right = force.cell_right
        trans_wall = calc_bc_wall_trans(model)

        T = state.WallTemperature[cell_right]
        T_bc = pars.T_a

        bc_src = -(trans_wall * (T - T_bc))
        acc[cell_right] -= bc_src
    end



end
