function initial_adsorbed_concentration(model; kwargs...)

    # Extract values from kwargs
    if haskey(kwargs, :Pressure) && haskey(kwargs, :y) && haskey(kwargs, :Temperature)
        p_init = kwargs[:Pressure]
        y_init = kwargs[:y]
        temperature_init = kwargs[:Temperature]
        R = model.system.p.R
        ncells = Jutul.number_of_cells(model.domain)

        # Ensure vectors have correct length for cell-wise operations
        p_vec = isa(p_init, Number) ? fill(p_init, ncells) : p_init
        temp_vec = isa(temperature_init, Number) ? fill(temperature_init, ncells) : temperature_init

        cTot = p_vec ./ (R * temp_vec)
        c = y_init' .* cTot

        q_init = map(1:ncells) do i
            qstar = compute_equilibrium(model.system, c[i,:], temp_vec[i])
        end
        q_init = stack(q_init) # Convert Vector of SVectors to Matrix

    else
        error("Need to set Pressure, y and Temperature in order to determine initial adsorbed concentration in the column")
    end
    return q_init
end

function setup_adsorption_model(system::AdsorptionSystem;
    ncells = 100
)
    constants = system.p
    dx = sqrt(pi*constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = mocca_domain(mesh, system)

    model = Jutul.SimulationModel(domain, system, general_ad = true)
    return model
end

function setup_adsorption_state(model; kwargs...)
    system = model.system
    g = JutulDarcy.physical_representation(model.data_domain)
    ncells = prod(g.dims)

    q_init = initial_adsorbed_concentration(model; kwargs...)

    state0 = Jutul.setup_state(model;
        AdsorbedConcentration = q_init,
        kwargs...)

    return state0
end

function setup_adsorption_parameters(model; kwargs...)
    system = model.system
    volumes = model.data_domain[:volumes]
    solid_volume = volumes * (1 - system.p.Φ)
    fluid_vol = volumes * system.p.Φ

    parameters = Jutul.setup_parameters(model;
        SolidVolume=solid_volume,
        FluidVolume=fluid_vol,
        kwargs...
    )

    return parameters
end


