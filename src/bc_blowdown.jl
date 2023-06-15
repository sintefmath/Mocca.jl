using Parameters

@with_kw struct BlowdownBC{T}
    PH::T
    PI::T
    λ::T
    cell_right::Int
    t_stage
end


function pressure_right(force::BlowdownBC, time)

    @show "blowdown" #DEBUG

    cycle_time = sum(force.t_stage)
    numcycles = max(0,time - cycle_time)
    step_end = cumsum(force.t_stage)
    t =  mod(time, cycle_time) - step_end[3]


    PH = force.PH
    PI = force.PI
    λ = force.λ

    return (PI + (PH - PI) * exp(-λ * t))
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::BlowdownBC,
    time,
)

    state = storage.state

    pars = model.system.p
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model)

    # right side
    begin
        cell_right = force.cell_right
        P = state.Pressure[cell_right]
        T = state.Temperature[cell_right]
        y = state.y[:, cell_right] 

        P_bc = pressure_right(force, time)

        q = -trans * mob * (P_bc - P)

        cTot = P / (T * R)

        for i in eachindex(y)
            c = y[i] * cTot
            mysource =  -(q * c)
            acc[i, cell_right] -= mysource
        end

    end
  
end


function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::BlowdownBC,
    time,
)

    state = storage.state


    pars = model.system.p
    ρ_g = pars.ρ_g
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model)

    # right side
    begin
        cell_right = force.cell_right
        P = state.Pressure[cell_right]

        C_pg = state.C_pg[cell_right]
        avm = state.avm[cell_right]        

        P_bc = pressure_right(force, time)


        q = -trans * mob * (P_bc - P)


        bc_src = -(q * P / R * C_pg * avm)
        acc[cell_right] -= bc_src

    end

end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::BlowdownBC,
    time,
)

    state = storage.state

    pars = model.system.p
 
    # right side
    begin
        cell_right = force.cell_right
        trans_wall = calc_bc_wall_trans(model)

        T = state.WallTemperature[cell_right]
        T_bc = pars.T_a

        bc_src = -(trans_wall * (T_bc - T))
        acc[cell_right] -= bc_src
    end
end
