function plot_outlet(case,result)

    # Get the substates and subtimesteps used inside the simulator for visualisation
    substates, subtimesteps = Jutul.expand_to_ministeps(result);

    # We plot primary variables at the outlet through time
    outlet_cell = size(result.states[1][:y][1,:],1)
    f_outlet = Mocca.plot_cell(substates, case.sim.model, subtimesteps, outlet_cell);

    return f_outlet
end
