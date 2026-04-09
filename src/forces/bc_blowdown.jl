using Parameters

"""
Blowdown boundary condition #TODO add more description!

# Fields:
* `PH`: High pressure [Pa]
* `PI`: Intermediate pressure [Pa]
* `λ`:  BC rampup parameter
* `stage_start`: absolute start time (computed in `setup_forces`)
"""
@with_kw struct BlowdownBC{T}
    PH::T
    PI::T
    λ::T
    stage_start::Float64 = 0.0
end

function pressure_right(force::BlowdownBC, time)
    t = time - force.stage_start
    return force.PI + (force.PH - force.PI) * exp(-force.λ * t)
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

    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_right = Jutul.number_of_cells(model.domain)
    trans = calc_bc_trans(model, state, cell_right)

    # right side
    begin
        P = state.Pressure[cell_right]
        T = state.Temperature[cell_right]
        y = state.y[:, cell_right] 

        P_bc = pressure_right(force, time)

        q = -trans * mob * (P_bc - P)

        cTot = P / (T * GAS_CONSTANT)

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

    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_right = Jutul.number_of_cells(model.domain)
    trans = calc_bc_trans(model, state, cell_right)

    # right side
    begin
        P = state.Pressure[cell_right]

        C_pg = state.C_pg[cell_right]
        avm = state.avm[cell_right]        

        P_bc = pressure_right(force, time)


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
    force::BlowdownBC,
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


function Jutul.vectorization_length(bc::BlowdownBC, variant)
    # PH::T
    # PI::T
    # λ::T
    return 3
end

function Jutul.vectorize_force!(v, bc::BlowdownBC, variant)
    if variant == :all
        names = [:PH, :PI, :λ]
        v[1] = bc.PH
        v[2] = bc.PI
        v[3] = bc.λ
    else
        error("Variant $variant not supported")
    end
    return (names = names, )
end

function Jutul.devectorize_force(bc::BlowdownBC, X, meta, variant)
    if variant == :all
        PH = X[1]
        PI = X[2]
        λ = X[3]
        return BlowdownBC(PH, PI, λ, bc.stage_start)
    else
        error("Variant $variant not supported")
    end
end
