using Test
using Mocca
using Jutul

# Breakthrough of the feed of the Mocca 0.1.0 constants, with no heat
# transfer to the wall, into a bed full of nitrogen at 1 bar, as a Process;
# returns the bed states
function run_dcb_simulation(; ncells = 200, t_ads = 5000.0, maxdt = 5000.0)
    p = haghpanah_breakthrough(haghpanah_constants(h_in = 0.0, h_out = 0.0); ncells = ncells, duration = t_ads, max_dt = maxdt)
    r = simulate(p; timestep_selector_cfg = (Bed = (y = 0.01, Temperature = 10.0, Pressure = 10.0),))
    return unit_states(r.states), p.flowsheet[:Bed], r.timesteps
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

    # From the Process (flowsheet) form of the breakthrough. The legacy
    # boundary conditions, whose references these replace, agreed to within
    # 8e-5 in every quantity
    ref = (
        y_outlet_CO2 = 0.1486563818654154,
        y_outlet_N2 = 0.851343618134585,
        T_outlet = 301.7194262101073,
        T_inlet = 298.14999483661705,
        T_max = 301.7194262101073,
        P_outlet = 100004.81235633386,
        P_inlet = 101884.30990920354,
        total_CO2_adsorbed = 774405.5802776606,
        breakthrough_time = 486.88859725807845,
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
