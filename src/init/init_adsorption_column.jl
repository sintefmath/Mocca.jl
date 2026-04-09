"""
    setup_process_model(system, domain::Jutul.DataDomain)

Build a `SimulationModel` from a prepared domain and system.
"""
setup_process_model(system, domain::Jutul.DataDomain) =
    Jutul.SimulationModel(domain, system, general_ad = true)

"""
    setup_process_model(system, constants::ConstantsStruct; ncells = 200, kwarg...)

Build a complete `SimulationModel` from a constants struct.
Computes permeability, dispersion, mesh and domain internally.
Keyword arguments are forwarded to `mocca_domain`.
"""
function setup_process_model(system, constants::ConstantsStruct; ncells = 200, kwarg...)
    porosity = constants.Φ
    perm = compute_permeability(porosity, constants.d_p)
    disp = compute_dispersion(constants.D_m, constants.V0_inter, constants.d_p)

    mesh = column_mesh(ncells, constants.L, constants.r_in)

    domain = mocca_domain(mesh;
        porosity = porosity, permeability = perm,
        r_in = constants.r_in, r_out = constants.r_out,
        diffusion_coefficient = disp, thermal_conductivity = constants.K_z,
        adsorbent_density = constants.ρ_s,
        adsorbent_heat_capacity = constants.C_ps,
        wall_density = constants.ρ_w,
        wall_heat_capacity = constants.C_pw,
        wall_conductivity = constants.K_w,
        fluid_viscosity = constants.fluid_viscosity,
        fluid_density = constants.ρ_g,
        inner_htc = constants.h_in,
        outer_htc = constants.h_out,
        ambient_temperature = constants.T_a,
        kwarg...)

    return setup_process_model(system, domain)
end

function initial_adsorbed_concentration(model, t_init, p_init, y_init)
    # Check that necessary state variables are provided
    !ismissing(t_init) || error("The key Temperature must be present to initialize AdsorbedConcentration")
    !ismissing(p_init) || error("The key Pressure must be present to initialize AdsorbedConcentration")
    !ismissing(y_init) || error("The key y must be present to initialize AdsorbedConcentration")
    ncells = Jutul.number_of_cells(model.domain)

    # Ensure vectors have correct length for cell-wise operations
    p_vec = isa(p_init, Number) ? fill(p_init, ncells) : p_init
    temp_vec = isa(t_init, Number) ? fill(t_init, ncells) : t_init

    cTot = p_vec ./ (GAS_CONSTANT * temp_vec)
    c = y_init' .* cTot

    q_init = map(1:ncells) do i
        qstar = compute_equilibrium(model.system.isotherm, c[i, :], temp_vec[i])
    end
    q_init = stack(q_init) # Convert Vector of SVectors to Matrix
    return q_init
end

function setup_process_state(model; kwargs...)
        # Check if a state is already provided, otherwise use kwargs
        vars = get(kwargs, :state0, kwargs)

        # If not provided, the initial adsorbed concentration is calculated from the other state variables
        if haskey(kwargs, :AdsorbedConcentration)
            q_init = vars[:AdsorbedConcentration]
        else
            T0 = get(vars, :Temperature, missing)
            p_init = get(vars, :Pressure, missing)
            y_init = get(vars, :y, missing)
            q_init = initial_adsorbed_concentration(model, T0, p_init, y_init)
        end

        state0 = Jutul.setup_state(model;
            AdsorbedConcentration = q_init,
            vars...)

    return state0
end

function setup_process_parameters(model; kwargs...)
    porosity = model.data_domain[:porosity]
    volumes = model.data_domain[:volumes]
    solid_volume = volumes .* (1 .- porosity)
    fluid_vol = volumes .* porosity

    parameters = Jutul.setup_parameters(model;
        SolidVolume=solid_volume,
        FluidVolume=fluid_vol,
        kwargs...
    )

    return parameters
end

function setup_forces(model, stage_times, bcs;
    num_cycles=1, max_dt = 1.0)

    cycle_time = sum(stage_times)
    step_end = cumsum(stage_times)

    timesteps = Float64[]
    sim_forces = []

    for j = 1:num_cycles
        for i in eachindex(stage_times)
            bc_i = bcs[i]
            if hasfield(typeof(bc_i), :stage_start)
                prev_end = i == 1 ? 0.0 : step_end[i-1]
                bc_i = reconstruct(bc_i, stage_start = (j - 1) * cycle_time + prev_end)
            end
            numsteps = stage_times[i] / max_dt
            append!(timesteps, repeat([max_dt], Int(floor(numsteps))))
            append!(sim_forces, repeat([Jutul.setup_forces(model, bc=bc_i)], Int(floor(numsteps))))
        end
    end

    return (sim_forces, timesteps)
end

"""
    setup_boundary_conditions(constants, stage_names)

Build boundary condition structs for each process stage from a constants struct.

`stage_names` is a vector of strings, one per stage. Recognised names:
`"pressurisation"`, `"adsorption"`, `"blowdown"`, `"evacuation"`.
"""
function setup_boundary_conditions(constants, stage_names)
    map(stage_names) do name
        if name == "pressurisation"
            PressurisationBC(
                y_feed = constants.y_feed, PH = constants.p_high, PL = constants.p_low,
                λ = constants.λ, T_feed = constants.T_feed)
        elseif name == "adsorption"
            AdsorptionBC(
                y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                T_feed = constants.T_feed)
        elseif name == "blowdown"
            BlowdownBC(
                PH = constants.p_high, PI = constants.p_intermediate,
                λ = constants.λ)
        elseif name == "evacuation"
            EvacuationBC(
                PL = constants.p_low, PI = constants.p_intermediate,
                λ = constants.λ)
        else
            error("Boundary condition type $name not recognized")
        end
    end
end
