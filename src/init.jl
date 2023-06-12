import MAT
"""
Read the initial data from a MATLAB setup. 

Return a simulator object.
"""
function initialize_from_matlab(datafilepath; general_ad::Bool=true, forcing_term_coefficient::Float64=1.0)
    data = MAT.matread(datafilepath)

    # TODO: Unify data files...
    if haskey(data, "problem")
        setup = data["problem"]["SimulatorSetup"]
        field = name -> data["problem"]["SimulatorSetup"]["state0"][name]
    elseif haskey(data, "problem_struct")
        setup = data["problem_struct"]["SimulatorSetup"]
        field = name -> data["problem_struct"]["SimulatorSetup"]["state0"][name]
    else
        setup = data["SimulatorSetup"]
        field = name -> data["SimulatorSetup"]["state0"][name]
    end

    parameters = read_adsorption_parameters_from_matlab(datafilepath)
    system = AdsorptionFlowSystem(forcing_term_coefficient=forcing_term_coefficient, p=parameters)
    perm = compute_permeability(system)
    flatten = x -> collect(Iterators.flatten(x))
    p = flatten(field("pressure"))
    temperature = flatten(field("T"))
    walltemperature = flatten(field("Twall"))
    qCO2 = field("qCO2")
    qN2 = field("qN2")
    q = hcat(qCO2, qN2)'
    yCO2 = field("yCO2")
    y = hcat(yCO2, 1.0 .- yCO2)'


    numberofcells = size(p, 1)

    @show p
    @info "Initializing simulator" p temperature walltemperature q y

    # This is major hack
    # TODO: How to get mesh length from matlab file?
    centroids = setup["model"]["G"]["cells"]["centroids"][:,1,1]
    approxdx = maximum(centroids[2:end]-centroids[1:end-1])
    extent = maximum(centroids .+ approxdx/2)

    # TODO: Remove this check
    @assert extent == 1.0

    dx = sqrt(pi*parameters.r_in^2)
    # TODO: Check the grid size. TODO: Don't hardcode extent.
    mesh = Jutul.CartesianMesh((numberofcells, 1, 1), (extent, dx, dx))

    domain = JutulDarcy.reservoir_domain(mesh, porosity=system.p.Φ, permeability=perm)
    domain[:diffusion_coefficient] = axial_dispersion(system)
    domain[:thermal_conductivity] = system.p.K_z # TODO: Check this and the one above

    model = Jutul.SimulationModel(domain, system, general_ad=general_ad)

    # TODO: Figure out a better way to compute the volumes
    volumes = ones(numberofcells) * prod(mesh.deltas)
    @assert volumes == domain[:volumes]
    solid_volume = volumes * (1 - system.p.Φ)
    fluid_volume = volumes * system.p.Φ

    parameters = Jutul.setup_parameters(model,
        solidVolume=solid_volume,
        fluidVolume=fluid_volume
    )

    state0 = Jutul.setup_state(model,
        Pressure=p,
        y=y,
        adsorptionRates=q,
        Temperature=temperature,
        WallTemperature=walltemperature)

    return (sim = Jutul.Simulator(model, state0=state0, parameters=parameters), state0 = state0, parameters = parameters)
end

function initialize_Haghpanah_model(datafilepath; general_ad::Bool=true, forcing_term_coefficient::Float64=1.0)


    parameters = haghpanah_parameters()
    system = AdsorptionFlowSystem(forcing_term_coefficient=forcing_term_coefficient, p=parameters)
    perm = compute_permeability(system)
    flatten = x -> collect(Iterators.flatten(x))
    p = flatten(field("pressure"))
    temperature = flatten(field("T"))
    walltemperature = flatten(field("Twall"))
    qCO2 = field("qCO2")
    qN2 = field("qN2")
    q = hcat(qCO2, qN2)'
    yCO2 = field("yCO2")
    y = hcat(yCO2, 1.0 .- yCO2)'


    numberofcells = size(p, 1)

    @info "Initializing simulator" p temperature walltemperature q y

    # This is major hack
    # TODO: How to get mesh length from matlab file?
    centroids = setup["model"]["G"]["cells"]["centroids"][:,1,1]
    approxdx = maximum(centroids[2:end]-centroids[1:end-1])
    extent = maximum(centroids .+ approxdx/2)

    # TODO: Remove this check
    @assert extent == 1.0

    dx = sqrt(pi*parameters.r_in^2)
    # TODO: Check the grid size. TODO: Don't hardcode extent.
    mesh = Jutul.CartesianMesh((numberofcells, 1, 1), (extent, dx, dx))

    domain = JutulDarcy.reservoir_domain(mesh, porosity=system.p.Φ, permeability=perm)
    domain[:diffusion_coefficient] = axial_dispersion(system)
    domain[:thermal_conductivity] = system.p.K_z # TODO: Check this and the one above

    model = Jutul.SimulationModel(domain, system, general_ad=general_ad)

    # TODO: Figure out a better way to compute the volumes
    volumes = ones(numberofcells) * prod(mesh.deltas)
    @assert volumes == domain[:volumes]
    solid_volume = volumes * (1 - system.p.Φ)
    fluid_volume = volumes * system.p.Φ

    parameters = Jutul.setup_parameters(model,
        solidVolume=solid_volume,
        fluidVolume=fluid_volume
    )

    state0 = Jutul.setup_state(model,
        Pressure=p,
        y=y,
        adsorptionRates=q,
        Temperature=temperature,
        WallTemperature=walltemperature)

    return (sim = Jutul.Simulator(model, state0=state0, parameters=parameters), state0 = state0, parameters = parameters)
end



