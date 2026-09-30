using Parameters

"""
Pressurisation boundary condition #TODO add more description!

# Fields:
* `y_feed`: Composition of feed gas
* `PH`: High pressure [Pa]
* `PL`: Low pressure [Pa]
* `λ`:  BC rampup parameter
* `T_feed`: Feed temperature [K]
* `stage_start`: absolute start time (computed in `setup_forces`)
"""
@with_kw struct PressurisationBC{T, N}
    y_feed::SVector{N,T}
    PH::T
    PL::T
    λ::T
    T_feed::T
    stage_start::Float64 = 0.0
end

function pressure_left(force::PressurisationBC, time)
    t = time - force.stage_start
    return force.PH - (force.PH - force.PL) * exp(-force.λ * t)
end

function mass_flux_left(state, model, time, force::PressurisationBC)
    mob = 1.0 / first(_state_or_domain(state, model, :FluidViscosity, :fluid_viscosity, Unit()))
    cell_left = 1
    trans = calc_bc_trans(model, state, cell_left)

    P = state.Pressure[cell_left]
    T = state.Temperature[cell_left]
    P_bc = pressure_left(force, time)
    y_bc = force.y_feed
    T_bc = force.T_feed

    # Volumetric flow out of the column, negative for inflow as in the other
    # boundary conditions
    q = -trans * mob * (P_bc - P)

    # Gas enters with the feed composition and leaves with that of the inlet
    # cell, so every component is conserved across the inlet (Danckwerts
    # condition) and no component goes negative when the flow reverses
    if q < 0
        c_tot = P_bc / (T_bc * GAS_CONSTANT)
        mass_flux = [q * y_bc[i] * c_tot for i in eachindex(y_bc)]
    else
        c_tot = P / (T * GAS_CONSTANT)
        mass_flux = [q * state.y[i, cell_left] * c_tot for i in eachindex(y_bc)]
    end
    return mass_flux
end

function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:ComponentMasses},
    eq_s,
    force::PressurisationBC,
    time,
)

    state = storage.state
    cell_left = 1

    mass_flux = mass_flux_left(state, model, time, force)
    for i in eachindex(mass_flux)
        acc[i, cell_left] += mass_flux[i]
    end

end


function Jutul.apply_forces_to_equation!(
    acc,
    storage,
    model::AdsorptionModel,
    eq::Jutul.ConservationLaw{:ColumnConservedEnergy},
    eq_s,
    force::PressurisationBC,
    time,
)

    state = storage.state

    ρ_g = state.FluidDensity[1]
    μ = state.FluidViscosity[1]
    mob = 1.0 / μ
    cell_left = 1
    trans = calc_bc_trans(model, state, cell_left)

    # left side
    begin
        P = state.Pressure[cell_left]
        T = state.Temperature[cell_left]
        C_pg = state.C_pg[cell_left]
        avm = state.AverageMolarMass[cell_left]

        P_bc = pressure_left(force, time)
        T_bc = force.T_feed

        q = -trans * mob * (P_bc - P)  # see mass_flux_left

        if q < 0
            # Inflow at the feed temperature, as in AdsorptionBC
            bc_src = -((q * ρ_g * C_pg * (T_bc - T)) + (q * P_bc / GAS_CONSTANT) * C_pg * avm)
        else
            # Outflow at the cell state, as in EvacuationBC
            bc_src = -(q * P / GAS_CONSTANT * C_pg * avm)
        end
        acc[cell_left] -= bc_src
    end
end


function Jutul.vectorization_length(bc::PressurisationBC, variant)
    # y_feed::SVector{N,T}
    # PH::T
    # PL::T
    # λ::T
    # T_feed::T
    return 4 + length(bc.y_feed)
end

function Jutul.vectorize_force!(v, bc::PressurisationBC, variant)
    if variant == :all
        names = [:PH, :PL, :λ, :T_feed]
        v[1] = bc.PH
        v[2] = bc.PL
        v[3] = bc.λ
        v[4] = bc.T_feed
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

function Jutul.devectorize_force(bc::PressurisationBC, X::AbstractVector{T}, meta, variant) where T
    if variant == :all
        PH = X[1]
        PL = X[2]
        λ = X[3]
        T_feed = X[4]
        N = length(bc.y_feed)
        tmp = zeros(T, N)
        for i = 1:N
            tmp[i] = X[i + 4]
        end
        y_feed = SVector{N, T}(tmp)
        return PressurisationBC(y_feed, PH, PL, λ, T_feed, bc.stage_start)
    else
        error("Variant $variant not supported")
    end
end
