function initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)
    system = model.system
    g = JutulDarcy.physical_representation(model.data_domain)
    ncells = prod(g.dims)
    R = model.system.p.R

    p_init = ones(ncells)*P_init
    temperature_init = ones(ncells)*T_init

    # assert(sum.(y_init)==1) #TODO figure out how to do this
    
    # TODO: check  how this works in mrst and change it
    cTot = p_init ./ (R * temperature_init)
    c = y_init .* cTot
    qN2 = ones(ncells)
    for i in 1:ncells
        qstar = compute_equilibrium(model.system, c[i,:], temperature_init[i])
        qN2[i] = qN2[i]*qstar[2]
    end
    qCO2 = ones(ncells)*0

    q_init = hcat(qCO2, qN2)

    walltemperature_init = ones(ncells)*Tw_init

    state0 = Jutul.setup_state(model,
        Pressure = p_init,
        y = y_init',
        AdsorbedConcentration = q_init',
        Temperature = temperature_init,
        WallTemperature = walltemperature_init)
    volumes = model.data_domain[:volumes]
    solid_volume = volumes * (1 - system.p.Φ)
    fluid_vol = volumes * system.p.Φ

    parameters = Jutul.setup_parameters(model,
        SolidVolume=solid_volume,
        FluidVolume=fluid_vol
    )
    return (state0, parameters)
end



