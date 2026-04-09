using Parameters


"""
Evacuation boundary condition #TODO add more description!

# Fields:
* `PL`: Low pressure [Pa]
* `PI`: Intermediate pressure [Pa]
* `λ`:  BC rampup parameter
* `stage_start`: absolute start time (computed in `setup_forces`)
"""
@with_kw struct EvacuationBC{T}
    PL::T
    PI::T
    λ::T
    stage_start::Float64 = 0.0
end

function pressure_left(force::EvacuationBC, time)
    t = time - force.stage_start
    return force.PL + (force.PI - force.PL) * exp(-force.λ * t)
end

function mass_flux_left(state, model, time, force::EvacuationBC)
    μ = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    cell_left = 1
    trans = calc_bc_trans(model, state, cell_left)

    P = state.Pressure[cell_left]
    T = state.Temperature[cell_left]
    y = state.y[:, cell_left]

    P_bc = pressure_left(force, time)
    q = -trans * mob * (P_bc - P)
    cTot = P / (T * GAS_CONSTANT)

    c = y .* cTot
    mass_flux = -q .* c
    return mass_flux
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
    cell_left = 1

    # left side
    mass_flux = mass_flux_left(state, model, time, force)
    for i in eachindex(mass_flux)
        acc[i, cell_left] -= mass_flux[i]
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

    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_left = 1
    trans = calc_bc_trans(model, state, cell_left)

    # left side
    begin
        P = state.Pressure[cell_left]

        C_pg = state.C_pg[cell_left]
        avm = state.avm[cell_left]

        P_bc = pressure_left(force, time)

        q = -trans * mob * (P_bc - P)

        bc_src = -(q * P / GAS_CONSTANT * C_pg * avm)
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

function Jutul.vectorization_length(bc::EvacuationBC, variant)
    # PL::T
    # PI::T
    # λ::T
    return 3
end

function Jutul.vectorize_force!(v, bc::EvacuationBC, variant)
    if variant == :all
        names = [:PL, :PI, :λ]
        v[1] = bc.PL
        v[2] = bc.PI
        v[3] = bc.λ
    else
        error("Variant $variant not supported")
    end
    return (names = names, )
end

function Jutul.devectorize_force(bc::EvacuationBC, X, meta, variant)
    if variant == :all
        PL = X[1]
        PI = X[2]
        λ = X[3]
        return EvacuationBC(PL, PI, λ, bc.stage_start)
    else
        error("Variant $variant not supported")
    end
end
