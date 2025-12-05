using Jutul

function plot_outlet(case, states, timesteps_out)


    outlet_cell = size(states[1][:y][1,:],1)
    f_outlet = Mocca.plot_cell(states, case.model, timesteps_out, outlet_cell);

    return f_outlet
end
