import MAT
"""
Read the initial data from a MATLAB setup. 

Return a simulator object.
"""
function initialize_from_matlab(datafilepath; general_ad::Bool=true, forcing_term_coefficient::Float64=1.0)
    data = MAT.matread(datafilepath)
    field = name -> data["problem"]["SimulatorSetup"]["state0"][name]

    # TODO: Read the relevant properties from the matlab file here (eg porosity)
    parameters = AdsorptionParameters()
    system = AdsorptionFlowSystem(forcing_term_coefficient=forcing_term_coefficient, p = parameters)
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

    nc = size(p, 1)
    p = collect(LinRange(p[1], 0.9*p[1], nc))

    @info "Initializing simulator" p temperature walltemperature q y

    G = JutulDarcy.get_1d_reservoir(nc, general_ad=general_ad, poro=system.p.Φ, perm=perm)

    model = Jutul.SimulationModel(G, system)

    g = Jutul.physical_representation(model.domain)

    pv = g.pore_volumes
    volumes = pv / system.p.Φ
    solid_volume = volumes * (1 - system.p.Φ)


    parameters = Jutul.setup_parameters(model,
        solidVolume=solid_volume,
    )

    state0 = Jutul.setup_state(model,
        Pressure=p,
        y=y,
        adsorptionRates=q,
        Temperature=temperature,
        WallTemperature=walltemperature)

    return Jutul.Simulator(model, state0=state0, parameters=parameters)
end




