"""
Adsorption boundary condition #TODO add more description!

# Fields:
* `y_feed`: Composition of feed gas
* `PH`: High pressure [Pa]
* `v_feed`: Feed velocity [m/s]
* `T_feed`: Feed temperature [K]
"""
@with_kw struct AdsorptionBC{T, N}
    y_feed::SVector{N,T}
    PH::T
    v_feed::T
    T_feed::T
end



function flux_left(model::AdsorptionModel, state, force::AdsorptionBC)
    Af = compute_column_face_area(model, state)
    return -force.v_feed * Af
end

function mass_flux_left(state, model, time, force::AdsorptionBC)
    μ = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    cell_left = 1
    trans = calc_bc_trans(model, state, cell_left)

    P = state.Pressure[cell_left]
    y = state.y[:, cell_left]
    y_bc = force.y_feed
    T_bc = force.T_feed

    q = flux_left(model, state, force)
    P_bc = q / (trans * mob) + P

    c_tot = P_bc / (T_bc * GAS_CONSTANT)
    c = y_bc .* c_tot

    mass_flux = c_tot .* q .* (y_bc .- y) .+ q .* c
    return mass_flux
end

function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    eq_s,
    force::AdsorptionBC,
    time,
)

    state = storage.state

    y = state.y
    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_left = 1
    cell_right = Jutul.number_of_cells(model.domain)
    trans = calc_bc_trans(model, state, cell_left)


    # left side
    mass_flux = mass_flux_left(state, model, time, force)
    for i in axes(y, 1)
        acc[i, cell_left] += mass_flux[i]
    end

    # right side
    @inbounds begin
        P = state.Pressure[cell_right]
        T = state.Temperature[cell_right]

        P_bc = force.PH
        y_bc = force.y_feed
        T_bc = force.T_feed

        q = -trans * mob * (P_bc - P)

        cTot = P / (T * GAS_CONSTANT)

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

    ρ_g = state.FluidDensity[1]
    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_left = 1
    cell_right = Jutul.number_of_cells(model.domain)
    trans = calc_bc_trans(model, state, cell_left)


    # left side
    @inbounds begin

        P = state.Pressure[cell_left]
        T = state.Temperature[cell_left]
        C_pg = state.C_pg[cell_left]
        avm = state.AverageMolarMass[cell_left]

        q = flux_left(model, state, force)

        P_bc = q / trans / mob + P
        y_bc = force.y_feed
        T_bc = force.T_feed


        cTot = P_bc / (T_bc * GAS_CONSTANT)
        c = y_bc .* cTot

        bc_src = -((q * ρ_g * C_pg * (T_bc - T)) + (q * P_bc / GAS_CONSTANT) * C_pg * avm)
        acc[cell_left] -= bc_src
    end

    # right side
    @inbounds begin

        P = state.Pressure[cell_right]
        C_pg = state.C_pg[cell_right]
        avm = state.AverageMolarMass[cell_right]

        P_bc = force.PH
        y_bc = force.y_feed
        T_bc = force.T_feed

        q = -trans * mob * (P_bc - P)

        bc_src = -(q * P / GAS_CONSTANT * C_pg * avm)
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

    T_bc = first(model.data_domain[:ambient_temperature, Column()])
    cell_left = 1
    cell_right = Jutul.number_of_cells(model.domain)

    # left side
    begin
        trans_wall = calc_bc_wall_trans(model, state, cell_left)

        T = state.WallTemperature[cell_left]

        bc_src = -(trans_wall * (T - T_bc))
        acc[cell_left] -= bc_src
    end

    # right side
    begin
        trans_wall = calc_bc_wall_trans(model, state, cell_right)

        T = state.WallTemperature[cell_right]

        bc_src = -(trans_wall * (T - T_bc))
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
        return AdsorptionBC(y_feed, PH, v_feed, T_feed)
    else
        error("Variant $variant not supported")
    end
end
