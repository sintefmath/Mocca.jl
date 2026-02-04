using Test
using Mocca
using Jutul

function run_dcb_simulation(; ncells = 200, t_ads = 5000.0, maxdt = 5000.0)
    constants = Mocca.HaghpanahConstants{Float64}(h_in = 0.0, h_out = 0.0)
    system = Mocca.TwoComponentAdsorptionSystem(constants)
    model = Mocca.setup_process_model(system; ncells = ncells)

    bar = si_unit(:bar)
    state0 = Mocca.setup_process_state(model;
        Pressure = 1 * bar,
        Temperature = 298.15,
        WallTemperature = constants.T_a,
        y = [1e-10, 1.0 - 1e-10]
    )

    parameters = Mocca.setup_process_parameters(model)
    sim_forces, timesteps = Mocca.setup_forces(model, [t_ads], ["adsorption"];
        num_cycles = 1, max_dt = maxdt)

    case = Mocca.MoccaCase(model, timesteps, sim_forces;
        state0 = state0, parameters = parameters)

    states, timesteps_out = Mocca.simulate_process(case;
        timestep_selector_cfg = (y = 0.01, Temperature = 10.0, Pressure = 10.0),
        output_substates = true,
        info_level = -1
    )

    return states, model, timesteps_out
end

function dcb_metrics(states, model, timesteps_out)
    final = states[end]
    nc = number_of_cells(model.domain)

    breakthrough_time = nothing
    t = 0.0
    for (i, state) in enumerate(states)
        t += timesteps_out[i]
        if state[:y][1, nc] > 0.01
            breakthrough_time = t
            break
        end
    end

    return (
        y_outlet_CO2 = final[:y][1, nc],
        y_outlet_N2 = final[:y][2, nc],
        T_outlet = final[:Temperature][nc],
        T_inlet = final[:Temperature][1],
        T_max = maximum(final[:Temperature]),
        P_outlet = final[:Pressure][nc],
        P_inlet = final[:Pressure][1],
        total_CO2_adsorbed = sum(final[:AdsorbedConcentration][1, :]),
        breakthrough_time = breakthrough_time,
    )
end

@testset "DCB Regression" begin
    states, model, timesteps_out = run_dcb_simulation()
    m = dcb_metrics(states, model, timesteps_out)

    ref = (
        y_outlet_CO2 = 0.1486736097381964,
        y_outlet_N2 = 0.851326390261803,
        T_outlet = 301.6723625454626,
        T_inlet = 298.1499948649097,
        T_max = 301.6723625454626,
        P_outlet = 100004.81147504161,
        P_inlet = 101884.19452854067,
        total_CO2_adsorbed = 774447.598037862,
        breakthrough_time = 484.5518204728325,
    )

    @test isapprox(m.y_outlet_CO2, ref.y_outlet_CO2, rtol=1e-3)
    @test isapprox(m.y_outlet_N2, ref.y_outlet_N2, rtol=1e-3)
    @test isapprox(m.T_outlet, ref.T_outlet, atol=0.5)
    @test isapprox(m.T_inlet, ref.T_inlet, atol=0.5)
    @test isapprox(m.T_max, ref.T_max, atol=0.5)
    @test isapprox(m.P_outlet, ref.P_outlet, rtol=1e-3)
    @test isapprox(m.P_inlet, ref.P_inlet, rtol=1e-3)
    @test isapprox(m.total_CO2_adsorbed, ref.total_CO2_adsorbed, rtol=1e-3)
    @test isapprox(m.breakthrough_time, ref.breakthrough_time, rtol=1e-3)
end
