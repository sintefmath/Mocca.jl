@with_kw struct AdsorptionBC{T, N}
    y_feed::SVector{N,T}
    PH::T
    v_feed::T
    T_feed::T
    cell_left::Int
    cell_right::Int
end



function flux_left(model::AdsorptionModel, force::AdsorptionBC)
    Af = compute_column_face_area(model)
    return -force.v_feed * Af
end


function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:TotalMasses},
    eq_s,
    force::AdsorptionBC,
    time,
)

    state = storage.state

    y = state.y
    pars = model.system.p
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model)


    # left side
    begin
        cell_left = force.cell_left
        P = state.Pressure[cell_left]

        q = flux_left(model,force)

        P_bc = q / trans / mob + P
        y_bc = force.y_feed        
        T_bc = force.T_feed

        cTot = P_bc / (T_bc * R)
        c = y_bc .* cTot

        for i in axes(y, 1)
            mysource = -(cTot * q * (y_bc[i] - y[i, cell_left]) + q * c[i])
            acc[i, cell_left] -= mysource
        end
    end

    # right side
    begin
        cell_right = force.cell_right
        P = state.Pressure[cell_right]
        T = state.Temperature[cell_right]

        P_bc = force.PH
        y_bc = force.y_feed
        T_bc = force.T_feed

        q = -trans * mob * (P_bc - P)

        cTot = P / (T * R)

        for i in axes(y, 1)
            c = y[i, cell_right] *cTot
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
    force::AdsorptionBC,
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
        T = state.Temperature[cell_left]
        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]

        q = flux_left(model,force)

        P_bc = q / trans / mob + P
        y_bc = force.y_feed
        T_bc = force.T_feed


        cTot = P_bc / (T_bc * pars.R)
        c = y_bc .* cTot

        bc_src = -((q * ρ_g * C_pg * (T_bc - T)) + (q * P_bc / R) * C_pg * avm)
        acc[cell_left] -= bc_src

    end
   
    # right side
    begin

        cell_right = force.cell_right
        P = state.Pressure[cell_right]
        C_pg = state.C_pg[cell_right]
        avm = state.avm[cell_right]

        P_bc = force.PH
        y_bc = force.y_feed
        T_bc = force.T_feed

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
    force::AdsorptionBC,
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

        bc_src = -(trans_wall * (T_bc - T))
        acc[cell_right] -= bc_src
    end

  
end