using Parameters

@with_kw struct EvacuationBC{T}
    PL::T
    PI::T
    λ::T
    cell_left::Int
end

pressure_function(PL, PI, λ, t) = (PI + (PL - PI) * exp(-λ * t))

function pressure_left(force::EvacuationBC, time)
    return pressure_function(force.PL, force.PI, force.λ, time)
end

function calc_bc_trans(model::AdsorptionFlowModel)
    k = Mocca.compute_permeability(model.system)
    dx = Mocca.compute_dx(model, 1) / 2
    A = (pi * model.system.p.r_in^2)
    return k * A / dx
end

function calc_bc_wall_trans(model::AdsorptionFlowModel)
    k = model.system.p.K_w
    dx = compute_dx(model, 1) / 2
    A = area_wall(model.system)
    return k * A / dx
end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
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

        q = -trans * mob * (P - P_bc)

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
    model::AdsorptionFlowModel,
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


        q = -trans * mob * (P - P_bc)


        bc_src = -(q * P / R * C_pg * avm)
        acc[cell_left] -= bc_src

    end

end




function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionFlowModel,
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
end
