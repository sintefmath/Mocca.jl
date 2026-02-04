using Test
using Mocca
using Jutul
using Statistics: mean

function run_cyclic_simulation(; ncells = 200, num_cycles = 3, max_dt = 1.0)
    constants = Mocca.HaghpanahConstants{Float64}()
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

    stage_times = [15.0, 15.0, 30.0, 40.0]
    stage_names = ["pressurisation", "adsorption", "blowdown", "evacuation"]
    sim_forces, timesteps = Mocca.setup_forces(model, stage_times, stage_names;
        num_cycles = num_cycles, max_dt = max_dt)

    case = Mocca.MoccaCase(model, timesteps, sim_forces;
        state0 = state0, parameters = parameters)

    states, timesteps_out = Mocca.simulate_process(case;
        output_substates = true,
        info_level = -1
    )

    return states, model, timesteps_out, num_cycles, stage_times
end

function cyclic_metrics(states, model, timesteps_out, num_cycles, stage_times)
    nc = number_of_cells(model.domain)
    final = states[end]

    # Find start of last cycle by tracking cumulative time
    cycle_duration = sum(stage_times)
    last_cycle_start_time = (num_cycles - 1) * cycle_duration

    # Find index where last cycle starts
    cumtime = 0.0
    last_cycle_start_idx = 1
    for i in 1:length(timesteps_out)
        if cumtime >= last_cycle_start_time
            last_cycle_start_idx = i
            break
        end
        cumtime += timesteps_out[i]
    end

    last_cycle_states = states[last_cycle_start_idx:end]

    return (
        y_CO2_final = final[:y][1, nc],
        y_N2_final = final[:y][2, nc],
        T_final = final[:Temperature][nc],
        P_final = final[:Pressure][nc],
        total_CO2_adsorbed_final = sum(final[:AdsorbedConcentration][1, :]),
        avg_T_last_cycle = mean([maximum(s[:Temperature]) for s in last_cycle_states]),
        avg_P_last_cycle = mean([mean(s[:Pressure]) for s in last_cycle_states]),
    )
end

@testset "Cyclic VSA Regression" begin
    states, model, timesteps_out, num_cycles, stage_times = run_cyclic_simulation()
    m = cyclic_metrics(states, model, timesteps_out, num_cycles, stage_times)

    ref = (
        y_CO2_final = 7.100524006314865e-8,
        y_N2_final = 0.9999999289947609,
        T_final = 295.76741982625623,
        P_final = 9998.130045421889,
        total_CO2_adsorbed_final = 52325.83782488843,
        avg_T_last_cycle = 348.8714798098076,
        avg_P_last_cycle = 40071.222347007766,
    )

    @test isapprox(m.y_CO2_final, ref.y_CO2_final, rtol=1e-3, atol=1e-10)
    @test isapprox(m.y_N2_final, ref.y_N2_final, rtol=1e-3)
    @test isapprox(m.T_final, ref.T_final, atol=0.5)
    @test isapprox(m.P_final, ref.P_final, rtol=1e-3)
    @test isapprox(m.total_CO2_adsorbed_final, ref.total_CO2_adsorbed_final, rtol=1e-3)
    @test isapprox(m.avg_T_last_cycle, ref.avg_T_last_cycle, atol=0.5)
    @test isapprox(m.avg_P_last_cycle, ref.avg_P_last_cycle, rtol=1e-3)
end
