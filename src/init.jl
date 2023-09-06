export initialise_state_AdsorptionColumn

function initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)

    ncells = model.ncells
    R = model.parameters.R

    p_init = ones(ncells)*P_init
    temperature_init = ones(ncells)*T_init

    assert(sum.(y_init)==1)
    
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

    walltemperature_init = ones(ncells)*parameters.T_a

    state0 = Jutul.setup_state(model,
        Pressure = p_init,
        y = y_init',
        AdsorbedConcentration = q_init',
        Temperature = temperature_init,
        WallTemperature = walltemperature_init)

    return state0
end



