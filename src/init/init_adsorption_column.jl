function initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)
    system = model.system
    g = JutulDarcy.physical_representation(model.data_domain)
    ncells = prod(g.dims)
    R = model.system.p.R

    p_init = fill(P_init, ncells)
    temperature_init = fill(T_init, ncells)

    cTot = p_init ./ (R * temperature_init)
    c = y_init .* cTot

    qN2 = map(1:ncells) do i
        qstar = compute_equilibrium(model.system, c[i,:], temperature_init[i])
        qstar[2]
    end

    qCO2 = zeros(eltype(qN2), ncells)

    q_init = hcat(qCO2, qN2)

    walltemperature_init = fill(Tw_init, ncells)

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



