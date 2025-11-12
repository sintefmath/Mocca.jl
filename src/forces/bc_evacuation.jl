using Parameters


"""
Evacuation boundary condition #TODO add more description!

# Fields:
* `PL`: Low pressure [Pa]
* `PI`: Intermediate pressure [Pa]
* `λ`:  BC rampup parameter
* `cell_left`: Cell where left bc is applied
* `cell_right`: Cell where right bc is applied
* `cycle_time`: Total time for one cycle
* `previous_step_end`: Time at the end of all previous cycle steps
"""

@with_kw struct EvacuationBC{T}
    PL::T
    PI::T
    λ::T
    cell_left::Int
    cell_right::Int
    cycle_time::Float64
    previous_step_end::Float64
end

function pressure_left(force::EvacuationBC, time)
    # NB: Since this is the last stage of the cycle,
    # some care is needed when using floor to determine the descrete cycle number.
    # Using a small ϵ to avoid that last time step of evacuation tips cycle_no.
    ϵ = 1e-6
    cycle_no = floor(time/(force.cycle_time+ϵ))

    t_0 = cycle_no*force.cycle_time + force.previous_step_end
    t = time - t_0

    PL = force.PL
    PI = force.PI
    λ = force.λ

    return (PL + (PI - PL) * exp(-λ * t))
end

function mass_flux_left(state, model, time, force::EvacuationBC)
    pars = model.system.p
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model, state)

    cell_left = force.cell_left
    P = state.Pressure[cell_left]
    T = state.Temperature[cell_left]
    y = state.y[:, cell_left]

    P_bc = pressure_left(force, time)
    q = -trans * mob * (P_bc - P)
    cTot = P / (T * R)

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
    cell_left = force.cell_left

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


    pars = model.system.p
    ρ_g = pars.ρ_g
    R = pars.R
    μ = pars.fluid_viscosity
    mob = 1.0 / μ
    trans = calc_bc_trans(model, state)

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
        return EvacuationBC(PL, PI, λ, bc.cell_left, bc.cell_right)
    else
        error("Variant $variant not supported")
    end
end
