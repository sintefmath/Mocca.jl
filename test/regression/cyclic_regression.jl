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

# The averages over the last cycle are weighted by time, so that they do not
# depend on where the solver splits a step
function cyclic_metrics(states, model, timesteps_out, num_cycles, stage_times)
    nc = number_of_cells(model.domain)
    final = states[end]
    w = Mocca.cycle_window(timesteps_out, sum(stage_times), num_cycles)
    time_average(f) = sum(f(states[i]) * timesteps_out[i] for i in w) / sum(timesteps_out[w])

    return (
        y_CO2_final = final[:y][1, nc],
        y_N2_final = final[:y][2, nc],
        T_final = final[:Temperature][nc],
        P_final = final[:Pressure][nc],
        total_CO2_adsorbed_final = sum(final[:AdsorbedConcentration][1, :]),
        avg_T_last_cycle = time_average(s -> maximum(s[:Temperature])),
        avg_P_last_cycle = time_average(s -> mean(s[:Pressure])),
    )
end

@testset "Cyclic VSA Regression" begin
    states, model, timesteps_out, num_cycles, stage_times = run_cyclic_simulation()
    m = cyclic_metrics(states, model, timesteps_out, num_cycles, stage_times)

    # From the Process (flowsheet) form of the cycle. The legacy boundary
    # conditions agreed with it to within 4e-4. Carrying the stream enthalpy
    # into the bed without an extra conduction term at the ports then lowered
    # the CO2 held by 1.1% and raised the average peak temperature by 0.5 K;
    # the scale-free convergence criteria moved nothing by more than 1e-4
    # apart from the trace of CO2 at the top of the bed. Computed at the
    # default nonlinear tolerance of 1e-6
    ref = (
        y_CO2_final = 6.564727010646069e-8,
        y_N2_final = 0.9999999343527307,
        T_final = 295.76761284065645,
        P_final = 9998.202956436902,
        total_CO2_adsorbed_final = 50242.140960161836,
        avg_T_last_cycle = 348.84602987208797,
        avg_P_last_cycle = 40166.56860054164,
    )

    @test isapprox(m.y_CO2_final, ref.y_CO2_final, rtol=1e-3, atol=1e-10)
    @test isapprox(m.y_N2_final, ref.y_N2_final, rtol=1e-3)
    @test isapprox(m.T_final, ref.T_final, atol=0.5)
    @test isapprox(m.P_final, ref.P_final, rtol=1e-3)
    @test isapprox(m.total_CO2_adsorbed_final, ref.total_CO2_adsorbed_final, rtol=1e-3)
    @test isapprox(m.avg_T_last_cycle, ref.avg_T_last_cycle, atol=0.5)
    @test isapprox(m.avg_P_last_cycle, ref.avg_P_last_cycle, rtol=1e-3)
end
