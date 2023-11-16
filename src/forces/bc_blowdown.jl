using Parameters

"""
Blowdown boundary condition #TODO add more description!

# Fields:
* `PH`: High pressure [Pa]
* `PI`: Intermediate pressure [Pa]
* `λ`:  BC rampup parameter
* `cell_left`: Cell where left bc is applied
* `cell_right`: Cell where right bc is applied
* `cycle_time`: Total time for one cycle
* `previous_step_end`: Time at the end of all previous cycle steps
"""

@with_kw struct BlowdownBC{T}
    PH::T
    PI::T
    λ::T
    cell_left::Int
    cell_right::Int
    cycle_time::Float64
    previous_step_end::Float64
end


function pressure_right(force::BlowdownBC, time)
    

    cycle_no = floor(time/force.cycle_time)

    t_0 = cycle_no*force.cycle_time + force.previous_step_end
    t = time - t_0

    PH = force.PH
    PI = force.PI
    λ = force.λ

    return (PI + (PH - PI) * exp(-λ * t))
end



function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
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
    trans = calc_bc_trans(model, state)

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
    model::AdsorptionModel,
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
    trans = calc_bc_trans(model, state)

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
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:WallConservedEnergy},
    eq_s,
    force::BlowdownBC,
    time,
)

    state = storage.state

    pars = model.system.p
 
    # left side
    begin
        cell_left = force.cell_left
        trans_wall = calc_bc_wall_trans(model, state)

        T = state.WallTemperature[cell_left]
        T_bc = pars.T_a

        bc_src = -(trans_wall * (T - T_bc))
        acc[cell_left] -= bc_src
    end

    # right side
    begin
        cell_right = force.cell_right
        trans_wall = calc_bc_wall_trans(model, state)

        T = state.WallTemperature[cell_right]
        T_bc = pars.T_a

        bc_src = -(trans_wall * (T - T_bc))
        acc[cell_right] -= bc_src
    end 
end
