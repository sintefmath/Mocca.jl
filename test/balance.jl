using Test
using Mocca
using Jutul

# Run a small DCB simulation for testing the balance functions
function _balance_test_sim(; ncells = 20, t_ads = 1000.0, maxdt = 100.0)
    constants = Mocca.HaghpanahConstants{Float64}(h_in = 0.0, h_out = 0.0)
    system = Mocca.AdsorptionSystem(constants)
    model = Mocca.setup_process_model(system, constants; ncells = ncells)

    bar = si_unit(:bar)
    state0 = Mocca.setup_process_state(model;
        Pressure      = 1 * bar,
        Temperature   = 298.15,
        WallTemperature = constants.T_a,
        y = [1e-10, 1.0 - 1e-10]
    )
    parameters = Mocca.setup_process_parameters(model)
    bcs = Mocca.setup_boundary_conditions(constants, ["adsorption"])
    sim_forces, timesteps = Mocca.setup_forces(model, [t_ads], bcs;
        num_cycles = 1, max_dt = maxdt)
    case = Mocca.MoccaCase(model, timesteps, sim_forces;
        state0 = state0, parameters = parameters)
    states, timesteps_out = Mocca.simulate_process(case;
        output_substates = true, info_level = -1)
    return (states, timesteps_out, model, sim_forces, state0)
end

@testset "Mass and Energy Balance" begin
    states, timesteps, model, forces, state0 = _balance_test_sim()

    @testset "mass_balance – full simulation" begin
        mb = Mocca.mass_balance(states, timesteps, model, forces; state0 = state0)
        # Should be small (good conservation); allow up to 5 % for coarse grid / large dt
        @test abs(mb) < 0.05
    end

    @testset "mass_balance – sub-period" begin
        t_half = sum(timesteps) / 2
        mb_half = Mocca.mass_balance(states, timesteps, model, forces;
            state0 = state0, t_end = t_half)
        @test abs(mb_half) < 0.05
    end

    @testset "mass_balance – component 2 (N₂)" begin
        mb_n2 = Mocca.mass_balance(states, timesteps, model, forces;
            state0 = state0, component = 2)
        @test abs(mb_n2) < 0.05
    end

    @testset "mass_balance – no state0" begin
        # Should still return a finite number without crashing
        mb = Mocca.mass_balance(states, timesteps, model, forces)
        @test isfinite(mb)
    end

    @testset "energy_balance – full simulation (adiabatic)" begin
        # h_in = h_out = 0 so no wall exchange; should conserve well
        eb = Mocca.energy_balance(states, timesteps, model, forces; state0 = state0)
        @test abs(eb) < 0.05
    end

    @testset "energy_balance – no state0" begin
        eb = Mocca.energy_balance(states, timesteps, model, forces)
        @test isfinite(eb)
    end

    @testset "energy_balance – sub-period" begin
        t_half = sum(timesteps) / 2
        eb_half = Mocca.energy_balance(states, timesteps, model, forces;
            state0 = state0, t_end = t_half)
        @test abs(eb_half) < 0.05
    end
end

@testset "Mass Balance – with wall heat transfer" begin
    # Non-adiabatic simulation (default h_in / h_out from constants)
    constants = Mocca.HaghpanahConstants{Float64}()   # h_in=8.6, h_out=2.5
    system = Mocca.AdsorptionSystem(constants)
    model = Mocca.setup_process_model(system, constants; ncells = 20)

    bar = si_unit(:bar)
    state0 = Mocca.setup_process_state(model;
        Pressure        = 1 * bar,
        Temperature     = 298.15,
        WallTemperature = constants.T_a,
        y = [1e-10, 1.0 - 1e-10]
    )
    parameters = Mocca.setup_process_parameters(model)
    bcs = Mocca.setup_boundary_conditions(constants, ["adsorption"])
    sim_forces, timesteps = Mocca.setup_forces(model, [500.0], bcs;
        num_cycles = 1, max_dt = 100.0)
    case = Mocca.MoccaCase(model, timesteps, sim_forces;
        state0 = state0, parameters = parameters)
    states, timesteps_out = Mocca.simulate_process(case;
        output_substates = true, info_level = -1)

    # Mass balance should still hold even with wall heat exchange
    mb = Mocca.mass_balance(states, timesteps_out, model, sim_forces; state0 = state0)
    @test abs(mb) < 0.05

    # Energy balance includes wall exchange; check that it returns a finite value
    eb = Mocca.energy_balance(states, timesteps_out, model, sim_forces; state0 = state0)
    @test isfinite(eb)
    @test abs(eb) < 0.10
end
