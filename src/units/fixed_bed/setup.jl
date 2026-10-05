"""
    column_mesh(ncells, L, r_in)

Create a 1D `CartesianMesh` for a cylindrical column.
The cross-section is a square with the same area as the circular column
(side length `√(π r_in²)`).
"""
function column_mesh(ncells, L, r_in)
    dx = sqrt(π * r_in^2)
    return Jutul.CartesianMesh((ncells, 1, 1), (L, dx, dx))
end

function compute_permeability(porosity, d_p)
    return 4 / 150 * ((porosity / (1 - porosity))^2) * (d_p / 2)^2 * porosity
end

function compute_dispersion(D_m, V0_inter, d_p)
    return 0.7 * D_m + 0.5 * V0_inter * d_p
end

"""
    mocca_domain(mesh; porosity, permeability, r_in, r_out, ...)

Set up a DataDomain for an adsorption column simulation.
Cell-level data (porosity, permeability) goes on Cells.
Column-level data (geometry, material, transport, heat transfer) goes on a
single Column entity.
"""
function mocca_domain(mesh;
    porosity, permeability,
    r_in, r_out,
    diffusion_coefficient, thermal_conductivity,
    adsorbent_density, adsorbent_heat_capacity,
    wall_density, wall_heat_capacity, wall_conductivity,
    fluid_viscosity,
    inner_htc, outer_htc, ambient_temperature,
    kwarg...
)
    domain = Jutul.DataDomain(mesh, porosity = porosity, permeability = permeability)

    # Register Column entity (count = 1)
    domain.entities[Unit()] = 1

    # Column-level geometry
    domain[:r_in, Unit()] = r_in
    domain[:r_out, Unit()] = r_out

    # Column-level transport coefficients
    domain[:diffusion_coefficient, Unit()] = diffusion_coefficient
    domain[:thermal_conductivity, Unit()] = thermal_conductivity

    # Column-level material properties
    domain[:adsorbent_density, Unit()] = adsorbent_density
    domain[:adsorbent_heat_capacity, Unit()] = adsorbent_heat_capacity
    domain[:wall_density, Unit()] = wall_density
    domain[:wall_heat_capacity, Unit()] = wall_heat_capacity
    domain[:wall_conductivity, Unit()] = wall_conductivity
    domain[:fluid_viscosity, Unit()] = fluid_viscosity
    domain[:inner_htc, Unit()] = inner_htc
    domain[:outer_htc, Unit()] = outer_htc
    domain[:ambient_temperature, Unit()] = ambient_temperature

    nc = Jutul.number_of_cells(mesh)
    dx = map(i -> first(Jutul.cell_dims(mesh, i)), 1:nc)
    domain[:dx] = dx

    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end

"""
    setup_process_model(system, domain::Jutul.DataDomain)

Build a `SimulationModel` from a prepared domain and system.
"""
setup_process_model(system, domain::Jutul.DataDomain) =
    Jutul.SimulationModel(domain, system, general_ad = true)

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

# A bed starts with its sorbent at equilibrium with the gas
function setup_unit_state(model::SorbentBedModel; pressure, temperature, composition, wall_temperature = temperature)
    return setup_process_state(model; _gas_state(model, pressure, temperature, composition, wall_temperature)...)
end

setup_unit_parameters(model::SorbentBedModel; kwarg...) = setup_process_parameters(model; kwarg...)
