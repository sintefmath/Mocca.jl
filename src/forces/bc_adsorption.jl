@with_kw struct AdsorptionBC{T, N}
    y_feed::SVector{N,T}
    PH::T
    v_feed::T
    T_feed::T
    cell_left::Int
    cell_right::Int
end



function flux_left(model::AdsorptionModel, state, force::AdsorptionBC)
    Af = compute_column_face_area(model, state)
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
    trans = calc_bc_trans(model, state)


    # left side
    @inbounds begin
        cell_left = force.cell_left
        P = state.Pressure[cell_left]

        q = flux_left(model, state, force)

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
    @inbounds begin
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
    trans = calc_bc_trans(model, state)


    # left side
    @inbounds begin

        cell_left = force.cell_left
        P = state.Pressure[cell_left]
        T = state.Temperature[cell_left]
        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]

        q = flux_left(model, state, force)

        P_bc = q / trans / mob + P
        y_bc = force.y_feed
        T_bc = force.T_feed


        cTot = P_bc / (T_bc * pars.R)
        c = y_bc .* cTot

        bc_src = -((q * ρ_g * C_pg * (T_bc - T)) + (q * P_bc / R) * C_pg * avm)
        acc[cell_left] -= bc_src

    end

    # right side
    @inbounds begin

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

        bc_src = -(trans_wall * (T_bc - T))
        acc[cell_right] -= bc_src
    end

  
end

function Jutul.vectorization_length(bc::AdsorptionBC, variant)
    # y_feed::SVector{N,T}
    # PH::T
    # v_feed::T
    # T_feed::T

    return 3 + length(bc.y_feed)
end

function Jutul.vectorize_force!(v, bc::AdsorptionBC, variant)
    if variant == :all
        names = [:PH, :v_feed, :T_feed]
        v[1] = bc.PH
        v[2] = bc.v_feed
        v[3] = bc.T_feed
        offset = length(names)
        for (i, f_i) in enumerate(bc.y_feed)
            offset += 1
            v[offset] = f_i
            push!(names, Symbol("y_feed$i"))
        end
    else
        error("Variant $variant not supported")
    end
    return (names = names, )
end

function Jutul.devectorize_force(bc::AdsorptionBC, X::AbstractVector{T}, meta, variant) where T
    if variant == :all
        PH = X[1]
        v_feed = X[2]
        T_feed = X[3]
        N = length(bc.y_feed)
        tmp = zeros(T, N)
        for i = 1:N
            tmp[i] = X[i + 3]
        end
        y_feed = SVector{N, T}(tmp)
        return AdsorptionBC(y_feed, PH, v_feed, T_feed, bc.cell_left, bc.cell_right)
    else
        error("Variant $variant not supported")
    end
end
