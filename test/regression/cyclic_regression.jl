using Test
using Mocca
using Jutul
using Statistics: mean

# The four-stage VSA of the Mocca 0.1.0 constants, from a bed full of nitrogen
# at 1 bar, as a Process; returns the bed states
function run_cyclic_simulation(; ncells = 200, num_cycles = 3, max_dt = 1.0)
    p = haghpanah_vsa(; ncells = ncells, max_dt = max_dt, initial_pressure = 1e5)
    r = simulate(p; cycles = num_cycles)
    stage_times = [st.duration for st in Mocca.stages(p)]
    return unit_states(r.states), p.flowsheet[:Bed], r.timesteps, num_cycles, stage_times
end

function cyclic_metrics(states, model, timesteps_out, num_cycles, stage_times)
    nc = number_of_cells(model.domain)
    final = states[end]

    cycle_duration = sum(stage_times)
    last_cycle_start_time = (num_cycles - 1) * cycle_duration

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

    # From the Process (flowsheet) form of the cycle. The legacy boundary
    # conditions, whose references these replace, agreed to within 4e-4 in
    # every quantity, apart from the trace of CO2 left at the top of the bed
    # (5.05e-8 there, 7.00e-8 here)
    ref = (
        y_CO2_final = 7.00076017689521e-8,
        y_N2_final = 0.9999999299923981,
        T_final = 295.76747069213474,
        P_final = 9998.232940271764,
        total_CO2_adsorbed_final = 50790.2436361158,
        avg_T_last_cycle = 348.2936828984303,
        avg_P_last_cycle = 40068.949813652536,
    )

    @test isapprox(m.y_CO2_final, ref.y_CO2_final, rtol=1e-3, atol=1e-10)
    @test isapprox(m.y_N2_final, ref.y_N2_final, rtol=1e-3)
    @test isapprox(m.T_final, ref.T_final, atol=0.5)
    @test isapprox(m.P_final, ref.P_final, rtol=1e-3)
    @test isapprox(m.total_CO2_adsorbed_final, ref.total_CO2_adsorbed_final, rtol=1e-3)
    @test isapprox(m.avg_T_last_cycle, ref.avg_T_last_cycle, atol=0.5)
    @test isapprox(m.avg_P_last_cycle, ref.avg_P_last_cycle, rtol=1e-3)
end
