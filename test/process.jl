using Test
using Mocca
using Jutul

@testset "Cycle metrics" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 20)
    mm = Mocca.setup_flowsheet_model(fs)
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.setup_process_state(fs[:Bed]; init...))
    prm = Mocca.setup_flowsheet_parameters(fs)
    forces, dt = Mocca.setup_schedule(fs, mm, stages; num_cycles = 2, max_dt = 1.0)
    states, ts = Mocca.simulate_process(Mocca.MoccaCase(mm, dt, forces; state0 = state0, parameters = prm);
        info_level = -1, output_substates = true)

    w2 = Mocca.cycle_window(ts, 100.0, 2)
    @test sum(ts[w2]) ≈ 100.0
    @test sum(ts[Mocca.cycle_window(ts, 100.0, 1)]) ≈ 100.0

    n_feed = Mocca.stream_totals(states, ts, :V_feed; window = w2)
    n_raff = Mocca.stream_totals(states, ts, :V_product; window = w2)
    n_ext = Mocca.stream_totals(states, ts, :V_vacuum; window = w2)
    # The feed stream has the feed composition
    @test n_feed[1] / sum(n_feed) ≈ constants.y_feed[1] rtol = 1e-6
    # Component balance over the cycle: feed − products = change in the bed
    inv_start = Mocca.component_inventory(fs[:Bed], states[first(w2) - 1][:Bed], prm[:Bed])
    inv_end = Mocca.component_inventory(fs[:Bed], states[last(w2)][:Bed], prm[:Bed])
    @test isapprox(n_feed .- n_raff .- n_ext, inv_end .- inv_start, rtol = 1e-6, atol = 1e-8 * sum(n_feed))

    p = Mocca.purity(n_ext)
    r = Mocca.recovery(n_ext, n_feed)
    @test 0 < p < 1
    @test 0 < r < 1
    @test p > constants.y_feed[1]  # the extract is enriched in CO2
    m_s = Mocca.sorbent_mass(fs[:Bed], prm[:Bed])
    @test m_s ≈ sum(prm[:Bed][:SolidVolume]) * constants.ρ_s
    @test Mocca.productivity(n_ext, 100.0, m_s) > 0

    W = Mocca.vacuum_pump_energy(states, ts, :V_vacuum; window = w2)
    @test W > 0
    e = Mocca.specific_energy_kwh_per_tonne(W, n_ext, constants.molecularMassOfCO2)
    @test isfinite(e) && e > 0
    # The feed valve never runs a vacuum pump
    @test Mocca.vacuum_pump_energy(states, ts, :V_feed; window = w2, discharge_pressure = 1.0) == 0
end

@testset "Cyclic steady state" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    init = (Pressure = constants.p_low, Temperature = 298.15, WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])

    # Legacy single-bed boundary conditions
    bed = Mocca.setup_process_model(Mocca.FixedBed(constants), constants; ncells = 10)
    state0 = Mocca.setup_process_state(bed; init...)
    prm = Mocca.setup_process_parameters(bed)
    bcs = Mocca.setup_boundary_conditions(constants, ["pressurisation", "adsorption", "blowdown", "evacuation"])
    forces, dt = Mocca.setup_forces(bed, [15.0, 15.0, 30.0, 40.0], bcs; num_cycles = 1, max_dt = 5.0)

    res = Mocca.simulate_to_cyclic_steady_state(bed, state0, prm, forces, dt; tol = 1e-3, max_cycles = 200)
    @test res.converged
    @test res.history[end] < 1e-3
    @test sum(res.timesteps) ≈ 100.0
    @test Mocca.cycle_change(bed, res.state0, Mocca.restart_state(bed, res.states[end])) ≈ res.history[end]

    # Successive substitution is much slower to get there
    plain = Mocca.simulate_to_cyclic_steady_state(bed, state0, prm, forces, dt; tol = 1e-3, max_cycles = res.cycles, acceleration = :none)
    @test !plain.converged
    @test plain.history[end] > res.history[end]

    # The same driver runs a flowsheet; at steady state the cycle's feed and
    # products balance
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 10)
    mm = Mocca.setup_flowsheet_model(fs)
    fs_state0 = Mocca.setup_flowsheet_state(fs; Bed = Mocca.setup_process_state(fs[:Bed]; init...))
    fs_prm = Mocca.setup_flowsheet_parameters(fs)
    fs_forces, fs_dt = Mocca.setup_schedule(fs, mm, stages; num_cycles = 1, max_dt = 5.0)
    fres = Mocca.simulate_to_cyclic_steady_state(mm, fs_state0, fs_prm, fs_forces, fs_dt; tol = 1e-3, max_cycles = 200)
    @test fres.converged
    n_feed = Mocca.stream_totals(fres.states, fres.timesteps, :V_feed)
    n_out = Mocca.stream_totals(fres.states, fres.timesteps, :V_product) .+ Mocca.stream_totals(fres.states, fres.timesteps, :V_vacuum)
    @test isapprox(sum(n_feed), sum(n_out), rtol = 1e-2)
end
