using Test
using Mocca
using Mocca.Units
using Jutul

const MU = Mocca

# A unit that is not a bed: an isothermal gas tank, a distributed unit with
# the gas block and nothing else. It only implements what the unit interface
# asks for (framework/unit_interface.jl), so a process holding it checks that the
# flowsheet and process layers do not assume every unit is a bed. `TS` says
# whether it has the TimeScale parameter the flowsheet requires.
struct _GasTank{TS} <: MU.DistributedUnit
    component_names::Vector{String}
    molecular_masses::Vector{Float64}
    thermodynamics::MU.IdealGasConstantCp{2, Float64}
    thermal::MU.Isothermal{Float64}
end

function Jutul.select_primary_variables!(S, ::_GasTank, model::Jutul.SimulationModel)
    S[:Pressure] = MU.Pressure()
    S[:y] = MU.GasMoleFractions()
    S[:Temperature] = MU.Temperature()
end

function Jutul.select_secondary_variables!(S, ::_GasTank, model::Jutul.SimulationModel)
    S[:TotalMolarConcentration] = MU.TotalMolarConcentration()
    S[:MolarConcentration] = MU.MolarConcentration()
    S[:AverageMolarMass] = MU.AverageMolarMass()
    S[:ComponentMasses] = MU.ComponentMasses()
end

function Jutul.select_equations!(eqs, sys::_GasTank, model::Jutul.SimulationModel)
    eqs[:mass_conservation] = Jutul.ConservationLaw(model.domain.discretizations.mass_flow, :ComponentMasses, MU.number_of_components(sys))
    eqs[:energy_column] = MU.FixedTemperatureEquation()
end

function Jutul.select_parameters!(S, ::_GasTank{TS}, model::Jutul.SimulationModel) where TS
    S[:FluidVolume] = MU.FluidVolume()
    S[:Permeability] = MU.Permeability()
    S[:CellDx] = MU.CellDx()
    S[:DispersionTransmissibilities] = MU.DispersionTransmissibilities()
    S[:AxialDispersion] = MU.AxialDispersion()
    S[:FluidViscosity] = MU.FluidViscosity()
    S[:BedCrossSectionArea] = MU.BedCrossSectionArea()
    S[:IsothermalTemperature] = MU.IsothermalTemperature()
    TS && (S[:TimeScale] = MU.TimeScale())
end

function _gas_tank_model(gas, temperature; volume, length = 1.0, time_scale = true)
    r = sqrt(volume / (π * length))
    mesh = MU.column_mesh(2, length, r)
    domain = Jutul.DataDomain(mesh; porosity = 1.0, permeability = 1e-6)
    domain.entities[MU.Unit()] = 1
    domain[:r_in, MU.Unit()] = r
    domain[:diffusion_coefficient, MU.Unit()] = 1e-5
    domain[:fluid_viscosity, MU.Unit()] = 1.8e-5
    domain[:dx] = fill(length / 2, 2)
    system = _GasTank{time_scale}(gas.names, gas.molar_masses, gas.thermodynamics, MU.Isothermal(temperature))
    return Jutul.SimulationModel(domain, system; general_ad = true)
end

# Moles of each component in the gas of a unit
_gas_inventory(state, parameters) = vec((state[:y] .* (state[:Pressure] ./ (MU.GAS_CONSTANT .* state[:Temperature]))') * parameters[:FluidVolume])

@testset "A unit that is not a bed in a process" begin
    gas = MU.haghpanah_2013_gas()
    feed = Feed(gas; composition = [0.15, 0.85], pressure = 1bar)
    fs = Flowsheet(gas, feed)
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), MU.haghpanah_2013_column(); ncells = 5)
    MU.add_unit!(fs, :Tank, _gas_tank_model(gas, feed.temperature; volume = 1.0))
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    # The bed's light product goes into the tank
    add_device!(fs, :V_link, :link; from = (:Bed, :top), to = (:Tank, :bottom))
    add_device!(fs, :V_product, :product; at = (:Tank, :top))
    @test MU.beds(fs) == [:Bed]
    @test MU.holdup_units(fs) == [:Bed, :Tank]
    @test MU.stream_equations(fs[:Tank]) == (:mass_conservation,)
    @test MU.stream_equations(fs[:Bed]) == (:mass_conservation, :energy_column)
    # Units dispatch to their own setup: a bed's solid volume excludes its pores
    prm = MU.setup_flowsheet_parameters(fs)
    @test prm[:Bed][:SolidVolume] == MU.setup_process_parameters(fs[:Bed])[:SolidVolume]
    @test prm[:Bed][:SolidVolume] ≈ fs[:Bed].data_domain[:volumes] .* (1 .- fs[:Bed].data_domain[:porosity])

    # Steps that name no bed act on the only bed
    p = Process(fs, Cycle([Pressurise(5.0), Adsorb(10.0; velocity = 0.5), Evacuate(10.0; to = 0.2bar)]);
        timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.1))
    @test isempty(validate(p))
    adsorption = MU.stages(p)[2].settings
    @test haskey(adsorption, :V_link)
    # The link's default conductance includes the tank's half cell
    C_bed, C_tank = MU.port_conductance(p.flowsheet[:Bed], 5), MU.port_conductance(p.flowsheet[:Tank], 1)
    @test adsorption[:V_link].law.conductance ≈ 1 / (1 / C_bed + 1 / C_tank)
    # Per connection, a port state cross term and a stream cross term for each
    # stream equation: three connections to the bed, two to the isothermal tank
    @test length(p.model.cross_terms) == 3 * (1 + 2) + 2 * (1 + 1)

    s0 = initial_state(p)
    @test all(s0[:Tank][:Pressure] .== 1bar) && !haskey(s0[:Tank], :AdsorbedConcentration)
    @test haskey(s0[:Bed], :AdsorbedConcentration)
    @test_throws MoccaInputError Process(fs, p.cycle; initial_state = Dict(:V_feed => (pressure = 1bar,)))

    # Every component is conserved over the bed, the tank and the streams
    r = simulate(p; state0 = s0)
    prm = p.parameters
    held(s) = MU.component_inventory(p.flowsheet[:Bed], s[:Bed], prm[:Bed]) .+ _gas_inventory(s[:Tank], prm[:Tank])
    n_in = MU.stream_totals(r.states, r.timesteps, :V_feed)
    n_out = MU.stream_totals(r.states, r.timesteps, :V_vacuum) .+ MU.stream_totals(r.states, r.timesteps, :V_product)
    @test sum(MU.stream_totals(r.states, r.timesteps, :V_link)) > 0
    @test isapprox(n_in .- n_out, held(r.states[end]) .- held(s0), rtol = 1e-6, atol = 1e-8 * sum(n_in))
    @test all(r.states[end][:Tank][:Temperature] .≈ feed.temperature)

    # Only the bed's sorbent counts in the performance
    @test performance(r).adsorbent_volume ≈ sum(prm[:Bed][:SolidVolume])
end

@testset "A unit without a time scale is refused" begin
    gas = MU.haghpanah_2013_gas()
    fs = MU.UnitNetwork()
    MU.add_unit!(fs, :Tank, _gas_tank_model(gas, 298.15; volume = 1.0, time_scale = false))
    MU.add_unit!(fs, :V, MU.setup_device_model(gas.names))
    MU.connect!(fs, (:Tank, :bottom) => (:V, :outlet))
    MU.set_boundary!(fs, :V, :inlet, MU.PortCondition(pressure = 1bar, temperature = 298.15, composition = [0.15, 0.85]))
    @test_throws ErrorException MU.setup_flowsheet_model(fs)
end
