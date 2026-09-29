@testset "System hierarchy and aliases" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.FixedBed(constants)

    @test Mocca.AdsorptionSystem === Mocca.FixedBed
    @test Mocca.Column === Mocca.Unit
    @test system isa Mocca.SorbentBed
    @test system isa Mocca.DistributedUnit
    @test system isa Mocca.MoccaSystem

    model = Mocca.setup_process_model(system, constants; ncells = 5)
    @test model isa Mocca.FixedBedModel
    @test model isa Mocca.MoccaModel
end

@testset "Boundary quantities are parameters" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.FixedBed(constants)
    model = Mocca.setup_process_model(system, constants; ncells = 5)
    parameters = Mocca.setup_process_parameters(model)

    @test haskey(parameters, :Permeability)
    @test haskey(parameters, :BedCrossSectionArea)
    @test parameters[:BedCrossSectionArea][1] ≈ π * constants.r_in^2
    @test all(parameters[:Permeability] .≈ Mocca.compute_permeability(constants.Φ, constants.d_p))

    # The half-cell transmissibility follows the parameter, not the data domain
    state = merge(Dict{Symbol, Any}(:CellDx => parameters[:CellDx]), parameters)
    state[:Permeability] = 2 .* parameters[:Permeability]
    nt = NamedTuple(state)
    @test Mocca.calc_bc_trans(model, nt, 1) ≈ 2 * Mocca.calc_bc_trans(model, parameters |> NamedTuple, 1)
end

function _short_breakthrough(system, constants; ncells = 20, t_end = 200.0, state_kw = (WallTemperature = constants.T_a,))
    model = Mocca.setup_process_model(system, constants; ncells = ncells)
    state0 = Mocca.setup_process_state(model;
        Pressure = constants.p_high,
        Temperature = constants.T_feed,
        y = [1e-10, 1.0 - 1e-10],
        state_kw...
    )
    parameters = Mocca.setup_process_parameters(model)
    bcs = Mocca.setup_boundary_conditions(constants, ["adsorption"])
    forces, dt = Mocca.setup_forces(model, [t_end], bcs; max_dt = 5.0)
    case = Mocca.MoccaCase(model, dt, forces; state0 = state0, parameters = parameters)
    states, = Mocca.simulate_process(case; info_level = -1)
    return model, states
end

@testset "Thermal models" begin
    constants = Mocca.HaghpanahConstants{Float64}()

    sys_wall = Mocca.FixedBed(constants)
    @test Mocca.has_wall(sys_wall)
    @test Mocca.has_energy_balance(sys_wall)

    sys_adiabatic = Mocca.FixedBed(constants; thermal = Mocca.Adiabatic())
    @test !Mocca.has_wall(sys_adiabatic)
    model_a, states_a = _short_breakthrough(sys_adiabatic, constants; state_kw = NamedTuple())
    @test !haskey(model_a.primary_variables, :WallTemperature)
    @test !haskey(model_a.equations, :energy_wall)
    @test maximum(states_a[end][:Temperature]) > constants.T_feed + 1.0

    sys_iso = Mocca.FixedBed(constants; thermal = Mocca.Isothermal(310.0))
    @test !Mocca.has_energy_balance(sys_iso)
    model_i, states_i = _short_breakthrough(sys_iso, constants; state_kw = NamedTuple())
    @test all(isapprox.(states_i[end][:Temperature], 310.0, atol = 1e-6))

    # Adiabatic bed runs hotter than one losing heat through the wall
    constants_hot = Mocca.HaghpanahConstants{Float64}(h_in = 50.0, h_out = 50.0)
    _, states_w = _short_breakthrough(Mocca.FixedBed(constants_hot), constants_hot)
    @test maximum(states_a[end][:Temperature]) > maximum(states_w[end][:Temperature])
end

@testset "Sorption heat over all components" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    sys_first = Mocca.FixedBed(constants)
    sys_all = Mocca.FixedBed(constants; sorption_heat_all_components = true)
    @test !Mocca.sorption_heat_all_components(sys_first)
    @test Mocca.sorption_heat_all_components(sys_all)

    # Check the term against a hand calculation for a two-component state
    q0 = [1.0 2.0; 3.0 4.0]
    q1 = [1.5 2.0; 2.0 4.0]  # cell 1: CO2 adsorbs, N2 desorbs
    st0 = (AdsorbedConcentration = q0,)
    st = (AdsorbedConcentration = q1, ΔH = [-30e3 -30e3; -15e3 -15e3], C_pa = [700.0, 700.0],
        AverageMolarMass = [0.03, 0.03], Temperature = [300.0, 300.0], SolidVolume = [0.1, 0.1])
    Δt = 2.0
    dq = (q1[:, 1] .- q0[:, 1]) ./ Δt
    per_component = [(700.0 * 0.03 * 300.0 + st.ΔH[i, 1]) * dq[i] for i in 1:2]
    @test Mocca.energy_source(Mocca.HeatOfSorption(), nothing, st, st0, 1, Δt) ≈ -0.1 * per_component[1]
    @test Mocca.energy_source(Mocca.HeatOfSorption(all_components = true), nothing, st, st0, 1, Δt) ≈ -0.1 * sum(per_component)
    @test Mocca.heat_capacity(Mocca.AdsorbedPhaseHeatCapacity(all_components = true), nothing, st, 1) ≈ 0.1 * 700.0 * 0.03 * (1.5 + 2.0)
    @test Mocca.heat_capacity(Mocca.AdsorbedPhaseHeatCapacity(), nothing, st, 1) ≈ 0.1 * 700.0 * 0.03 * 1.5

    # And the option reaches the simulation
    _, s_first = _short_breakthrough(sys_first, constants)
    _, s_all = _short_breakthrough(sys_all, constants)
    @test !isapprox(maximum(s_all[end][:Temperature]), maximum(s_first[end][:Temperature]), atol = 0.1)
end

# A unit can be extended by adding a source term, without touching the residuals
struct _ConstantHeater <: Mocca.AbstractSourceTerm
    power::Float64
end
Mocca.energy_source(h::_ConstantHeater, model, state, state0, cell, Δt) = h.power

@testset "Custom source term" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.FixedBed(constants; thermal = Mocca.Adiabatic())
    model = Mocca.setup_process_model(system, constants; ncells = 4)
    state0 = Mocca.setup_process_state(model; Pressure = 1e5, Temperature = 298.15, y = [1e-10, 1.0 - 1e-10])
    parameters = Mocca.setup_process_parameters(model)

    # No flow and no uptake: only the heater changes the temperature
    terms = (Mocca.SolidHeatCapacity(), Mocca.AdsorbedPhaseHeatCapacity())
    local_state = merge(NamedTuple(state0), NamedTuple(parameters))
    C = sum(Mocca.heat_capacity(t, model, local_state, 1) for t in terms)
    @test C > 0
    @test Mocca.sum_terms(Mocca.energy_source, (_ConstantHeater(5.0), _ConstantHeater(2.0)), model, local_state, local_state, 1, 1.0) == 7.0
    @test Mocca.sum_terms(Mocca.energy_source, (), model, local_state, local_state, 1, 1.0) == 0
end
