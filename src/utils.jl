export mocca_domain
"""
    mocca_domain(mesh::Jutul.CartesianMesh, system::AdsorptionSystem)

Set up a data domain for carbon capture simulation.
"""
function mocca_domain(mesh, system::AdsorptionSystem; kwarg...)
    domain = JutulDarcy.reservoir_domain(mesh, porosity = system.p.Φ, permeability = system.permeability)
    domain[:diffusion_coefficient] = system.dispersion
    domain[:thermal_conductivity] = system.p.K_z
    nc = Jutul.number_of_cells(mesh)

    dx = map(i -> first(Jutul.cell_dims(mesh, i)), 1:nc)
    domain[:dx] = dx

    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end

# TODO: This causes problem for adjoint simulation. Is it needed?
function Jutul.select_linear_solver(model::AdsorptionModel; kwarg...)
    #return Jutul.LUSolver(; kwarg...)
    return nothing
end

function simulate_adsorption(model, state0, dt, parameters, sim_forces;
    info_level = 0,
    var_tstep_cfg = nothing
)
    extra_args = Dict()

    # Set up optional variable timestep selectors
    if !isnothing(var_tstep_cfg)
        t_base = Jutul.TimestepSelector(initial_absolute = 1.0)
        timesteppers = Any[t_base]
        push!(timesteppers, t_base)
        for (k, v) in pairs(var_tstep_cfg)
            t_i = Jutul.VariableChangeTimestepSelector(k, v, relative = false)
            push!(timesteppers, t_i)
        end
        extra_args[:timestep_selectors] = timesteppers
    end

    sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)

    cfg = Jutul.simulator_config(sim;
        output_substates = true,
        info_level = info_level,
        pairs(extra_args)...
    )

    result = Jutul.simulate!(sim, dt;
        config = cfg,
        forces = sim_forces
    )

    substates, subtimesteps = Jutul.expand_to_ministeps(result)
    return substates, subtimesteps
end
