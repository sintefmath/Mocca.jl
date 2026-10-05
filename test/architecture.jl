@testset "System hierarchy and aliases" begin
    model = haghpanah_bed(ncells = 5)
    system = model.system

    @test Mocca.AdsorptionSystem === Mocca.FixedBed
    @test system isa Mocca.SorbentBed
    @test system isa Mocca.DistributedUnit
    @test system isa Mocca.MoccaSystem
    @test model isa Mocca.FixedBedModel
    @test model isa Mocca.MoccaModel
end

@testset "Boundary quantities are parameters" begin
    constants = haghpanah_constants()
    model = haghpanah_bed(constants; ncells = 5)
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

# A 200 s breakthrough with steps of 5 s; returns the bed model and its states
function _short_breakthrough(constants = haghpanah_constants(); ncells = 20, t_end = 200.0, kwarg...)
    p = haghpanah_breakthrough(constants; ncells = ncells, duration = t_end, max_dt = 5.0, kwarg...)
    r = simulate(p)
    return p.flowsheet[:Bed], unit_states(r.states)
end

@testset "Thermal models" begin
    constants = haghpanah_constants()

    sys_wall = haghpanah_bed(ncells = 5).system
    @test Mocca.has_wall(sys_wall)
    @test Mocca.has_energy_balance(sys_wall)

    model_a, states_a = _short_breakthrough(; thermal = Mocca.Adiabatic())
    @test !Mocca.has_wall(model_a.system)
    @test !haskey(model_a.primary_variables, :WallTemperature)
    @test !haskey(model_a.equations, :energy_wall)
    @test maximum(states_a[end][:Temperature]) > constants.T_feed + 1.0

    model_i, states_i = _short_breakthrough(; thermal = Mocca.Isothermal(310.0))
    @test !Mocca.has_energy_balance(model_i.system)
    @test all(isapprox.(states_i[end][:Temperature], 310.0, atol = 1e-6))

    # Adiabatic bed runs hotter than one losing heat through the wall
    _, states_w = _short_breakthrough(haghpanah_constants(h_in = 50.0, h_out = 50.0))
    @test maximum(states_a[end][:Temperature]) > maximum(states_w[end][:Temperature])
end

@testset "Sorption heat over all components" begin
    sys_first = haghpanah_bed(ncells = 5).system
    sys_all = haghpanah_bed(ncells = 5; sorption_heat = :all).system
    @test !Mocca.sorption_heat_all_components(sys_first)
    @test Mocca.sorption_heat_all_components(sys_all)

    # Check the term against a hand calculation for a two-component state
    q0 = [1.0 2.0; 3.0 4.0]
    q1 = [1.5 2.0; 2.0 4.0]  # cell 1: CO2 adsorbs, N2 desorbs
    st0 = (AdsorbedLoading = q0,)
    st = (AdsorbedConcentration = q1, AdsorbedLoading = q1, ΔH = [-30e3 -30e3; -15e3 -15e3], C_pa = [700.0, 700.0],
        AverageMolarMass = [0.03, 0.03], Temperature = [300.0, 300.0], SolidVolume = [0.1, 0.1])
    Δt = 2.0
    dq = (q1[:, 1] .- q0[:, 1]) ./ Δt
    per_component = [(700.0 * 0.03 * 300.0 + st.ΔH[i, 1]) * dq[i] for i in 1:2]
    @test Mocca.energy_source(Mocca.HeatOfSorption(), nothing, st, st0, 1, Δt) ≈ -0.1 * per_component[1]
    @test Mocca.energy_source(Mocca.HeatOfSorption(all_components = true), nothing, st, st0, 1, Δt) ≈ -0.1 * sum(per_component)
    @test Mocca.heat_capacity(Mocca.AdsorbedPhaseHeatCapacity(all_components = true), nothing, st, 1) ≈ 0.1 * 700.0 * 0.03 * (1.5 + 2.0)
    @test Mocca.heat_capacity(Mocca.AdsorbedPhaseHeatCapacity(), nothing, st, 1) ≈ 0.1 * 700.0 * 0.03 * 1.5

    # And the option reaches the simulation
    _, s_first = _short_breakthrough()
    _, s_all = _short_breakthrough(; sorption_heat = :all)
    @test !isapprox(maximum(s_all[end][:Temperature]), maximum(s_first[end][:Temperature]), atol = 0.1)
end

# A unit can be extended by adding a source term, without touching the residuals
struct _ConstantHeater <: Mocca.AbstractSourceTerm
    power::Float64
end
Mocca.energy_source(h::_ConstantHeater, model, state, state0, cell, Δt) = h.power

@testset "Custom source term" begin
    model = haghpanah_bed(ncells = 4, thermal = Mocca.Adiabatic())
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

@testset "Convergence criteria do not depend on the size of a cell" begin
    model = haghpanah_bed(ncells = 2)
    crit(eq, state, r; dt = 2.0) = Jutul.convergence_criterion(model, (state = state,), model.equations[eq], nothing, r; dt = dt)
    gas = (TotalMolarConcentration = [40.0, 20.0], FluidVolume = [1e-3, 2e-3])
    r = [1e-3 2e-3; 0.0 4e-4]
    # Mass: the residual over a step as a fraction of the gas in the cell
    c = crit(:mass_conservation, gas, r)
    @test c.CNV.errors ≈ [2e-3 * 2 / 0.04, 4e-4 * 2 / 0.04]
    @test c.CNV.names == ["CNV_CO2", "CNV_N2"]
    # A cell ten times larger with ten times the residual converges alike
    big = crit(:mass_conservation, merge(gas, (FluidVolume = 10 .* gas.FluidVolume,)), 10 .* r)
    @test big.CNV.errors ≈ c.CNV.errors

    # Energy: the temperature error that the residual over a step amounts to
    state = merge(gas, (C_pg = [1000.0, 1000.0], AverageMolarMass = [0.03, 0.03], SolidVolume = [2e-3, 4e-3],
        AdsorbentDensity = [1000.0], AdsorbentHeatCapacity = [900.0], AdsorbedConcentration = zeros(2, 2), C_pa = [800.0, 800.0]))
    C = [2e-3 * 1000 * 900 + 0.04 * 1000 * 0.03, 4e-3 * 1000 * 900 + 0.04 * 1000 * 0.03]
    e = crit(:energy_column, state, [5.0 -20.0])
    @test e.ΔT.errors ≈ [max(5.0 * 2 / C[1], 20.0 * 2 / C[2])]
end

@testset "Device equations converge relative to their flow" begin
    dev = Mocca.setup_device_model(["CO2", "N2"])
    tol = Jutul.set_default_tolerances(dev; tol = 1e-4)
    @test tol[:default] == 1e-4
    @test all(tol[k] ≈ 1e-7 for k in (:inlet_state, :outlet_state, :device_law))
    R = Mocca.GAS_CONSTANT
    st = (InletPressure = [2e5], OutletPressure = [1e5], InletTemperature = [300.0], OutletTemperature = [350.0], MolarFlow = [0.0])
    # A valve's flow scale is its flow with the whole upstream pressure as the
    # drop, so it stays finite when the pressures meet
    v = Mocca.LinearValve(1e-6)
    @test Mocca.flow_scale(v, st, 0.0) ≈ 1e-6 * 2e5^2 / (R * 300.0)
    @test Mocca.flow_scale(v, merge(st, (OutletPressure = [2e5],)), 0.0) ≈ 1e-6 * 2e5^2 / (R * 300.0)
    @test Mocca.flow_scale(Mocca.VolumetricFlow(0.01), st, 0.0) ≈ Mocca.device_flow(Mocca.VolumetricFlow(0.01), st, 0.0)
    # A device just closed measures what is left of its flow against itself
    @test Mocca.flow_scale(Mocca.Closed(), st, 0.0) == Mocca.MIN_FLOW_SCALE
    @test Mocca.flow_scale(Mocca.Closed(), merge(st, (MolarFlow = [-0.3],)), 0.0) == 0.3
end
