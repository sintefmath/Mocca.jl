using Test
using Mocca
using Mocca.Units
using Jutul

const M = Mocca

# The largest difference between the vectorised forces of two schedules
function force_difference(fa, ma, fb, mb)
    length(fa) == length(fb) || return Inf
    return maximum(maximum(abs.(first(Jutul.vectorize_forces(a, ma)) .- first(Jutul.vectorize_forces(b, mb))); init = 0.0) for (a, b) in zip(fa, fb))
end

# A step type defined outside Mocca: open the feed valve and hold the bed at
# the feed pressure
Base.@kwdef mutable struct HoldAtFeed <: Mocca.AbstractStep
    duration::Float64
    bed::Union{Nothing, Symbol} = nothing
    name::String = ""
end
Mocca.default_step_name(::HoldAtFeed) = "hold"
function Mocca.step_settings(s::HoldAtFeed, ctx::Mocca.StepContext)
    d = Mocca.find_device(ctx, (:feed,))
    return [d => (law = Mocca.valve(ctx, d), inlet = Mocca.boundary(ctx, d, ctx.fs.feed.pressure))]
end
Mocca.end_pressure(::HoldAtFeed, ctx) = ctx.fs.feed.pressure
register_step!("hold_at_feed", HoldAtFeed)

@testset "Units" begin
    @test 0.1bar == 1e4
    @test -36.6kJ_per_mol == -36600.0
    @test celsius(25) ≈ 298.15
    @test M.Units.parse_quantity("0.1 bar") == 1e4
    @test M.Units.parse_quantity("25 °C") ≈ 298.15
    @test M.Units.parse_quantity("100 cm/s") ≈ 1.0
    @test M.Units.parse_quantity("1e-3") == 1e-3
    @test_throws ArgumentError M.Units.parse_quantity("1 furlong")
end

@testset "Gas, feed and compositions" begin
    gas = Gas(["CO2", "N2", "H2O"]; heat_capacities = (H2O = 2000.0,))
    @test gas.heat_capacities == [846.0, 1040.0, 2000.0]
    @test composition(gas, (CO2 = 0.15, N2 = 0.85)) ≈ [0.15 * (1 - 1e-10), 0.85 * (1 - 1e-10), 1e-10]
    @test composition(gas, [0.1, 0.8, 0.1]) == [0.1, 0.8, 0.1]
    @test_throws MoccaInputError composition(gas, (CO2 = 0.15, N2 = 0.8))
    @test_throws MoccaInputError composition(gas, (CO2 = 0.15, Ne = 0.85))
    @test_throws MoccaInputError Gas(["CO2", "Xe"])
    register_component!(Component("Xe"; molar_mass = 131.29e-3, heat_capacity = 158.0, viscosity = 2.3e-5))
    @test Gas(["CO2", "Xe"]).molar_masses[2] == 131.29e-3
    feed = Feed(gas; composition = (CO2 = 0.15, N2 = 0.82, H2O = 0.03), pressure = 1bar)
    @test feed.light_product[1] == 1e-10
    @test sum(feed.light_product) ≈ 1
    # Viscosity and density default to those of the feed
    @test 1.4e-5 < M.gas_viscosity(gas, feed.composition) < 1.8e-5
    @test M.gas_density(gas, feed.composition, 1e5, 298.15) ≈ 1e5 * sum(feed.composition .* gas.molar_masses) / (M.GAS_CONSTANT * 298.15)
end

@testset "Isotherms and adsorbents" begin
    c = M.HaghpanahConstants{Float64}()
    old = DualSiteLangmuir(qsb = c.qsbi, b0 = c.b0, ΔUb = c.ΔUbi, qsd = c.qsdi, d0 = c.d0, ΔUd = c.ΔUdi)
    a = load_adsorbent("zeolite_13x_haghpanah_2013")
    gas = M.haghpanah_2013_gas()
    iso = M.isotherm(a, gas)
    @test iso.qsb == old.qsb && iso.b0 == old.b0 && iso.ΔUb == old.ΔUb && iso.ΔH == old.ΔH
    # Order follows the gas; components without data do not adsorb
    rev = M.isotherm(a, Gas(["N2", "H2O", "CO2"]))
    @test rev.qsb == [old.qsb[2], 0.0, old.qsb[1]]
    q = equilibrium_loading(a, gas; pressure = 1bar, temperature = 298.15, composition = (CO2 = 0.15, N2 = 0.85))
    @test 2 < q[1] < 4 && 0 < q[2] < q[1]
    tab = isotherm_table(a, gas, "CO2"; pressures = [0.01bar, 0.1bar, 1bar], temperatures = [298.15, 348.15])
    @test size(tab) == (3, 2) && issorted(tab[:, 1]) && all(tab[:, 2] .< tab[:, 1])
    # Capacities per kg are converted with the density
    perkg = dual_site_langmuir("CO2" => (qsb = 3.0, b0 = 1e-6, ΔUb = -30kJ_per_mol); basis = :mass)
    @test M.build_isotherm(perkg, Gas(["CO2", "N2"]); density = 1000.0).qsb == [3000.0, 0.0]
    @test langmuir("CO2" => (qs = 3.0, b0 = 1e-6, ΔU = -30kJ_per_mol)).qsd == [0.0]
    @test_throws MoccaInputError dual_site_langmuir("CO2" => (qsb = 1.0, b0 = 1e-6))
    @test_throws MoccaInputError dual_site_langmuir("CO2" => (qsb = 1.0, b0 = 1e-6, ΔUb = -1.0, qs = 2.0))
    @test_throws MoccaInputError Adsorbent(; isotherm = perkg, density = 1000.0, heat_capacity = 900.0, particle_diameter = 1e-3, particle_porosity = 1.5)
    @test_throws MoccaInputError Column(; length = 1.0, radius = 0.1, porosity = 0.0)
    @test library_adsorbents() == ["silica_gel_krishnamurthy_2014", "zeolite_13x_haghpanah_2013", "zeolite_13x_krishnamurthy_2014"]
    # Writing and reading an adsorbent gives it back
    path = save_adsorbent(joinpath(mktempdir(), "a.json"), a)
    b = load_adsorbent(path)
    @test M.isotherm(b, gas).ΔH == iso.ΔH && b.adsorbed_heat_capacity == a.adsorbed_heat_capacity
end

@testset "The old constants convert to the same bed" begin
    c = M.HaghpanahConstants{Float64}()
    k = convert_constants(c)
    fs = Flowsheet(k.gas, k.feed)
    add_bed!(fs, :Bed, k.adsorbent, k.column; ncells = 10, dispersion_velocity = k.dispersion_velocity)
    bed = fs[:Bed]
    iso = DualSiteLangmuir(qsb = c.qsbi, b0 = c.b0, ΔUb = c.ΔUbi, qsd = c.qsdi, d0 = c.d0, ΔUd = c.ΔUdi)
    @test bed.system.isotherm.qsb == iso.qsb && bed.system.isotherm.ΔH == iso.ΔH
    d = bed.data_domain
    @test all(d[:porosity] .== c.Φ) && all(d[:permeability] .== M.compute_permeability(c.Φ, c.d_p))
    @test first(d[:diffusion_coefficient, M.Unit()]) == M.compute_dispersion(c.D_m, c.V0_inter, c.d_p)
    @test first(d[:r_in, M.Unit()]) == c.r_in && first(d[:r_out, M.Unit()]) == c.r_out
    @test first(d[:adsorbent_density, M.Unit()]) == c.ρ_s && first(d[:wall_heat_capacity, M.Unit()]) == c.C_pw
    @test k.velocity ≈ c.v_feed / c.Φ && k.p_low == c.p_low
end

@testset "A component the sorbent does not take up" begin
    # The adsorbent has no O2 data, so O2 does not adsorb, and the bed
    # simulates with its loading staying at zero
    gas = Gas(["CO2", "N2", "O2"])
    fs = Flowsheet(gas, Feed(gas; composition = (CO2 = 0.15, N2 = 0.8, O2 = 0.05), pressure = 1bar))
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column(); ncells = 5)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    p = Process(fs, Cycle([Pressurise(5.0), Adsorb(10.0; superficial_velocity = 0.37), Evacuate(10.0; to = 0.1bar)]))
    r = simulate(p)
    bed = r.states[end][:Bed]
    @test all(isfinite, bed[:AdsorbedConcentration]) && all(isfinite, bed[:Pressure])
    @test maximum(bed[:AdsorbedConcentration][3, :]) < 1e-8
    @test maximum(bed[:AdsorbedConcentration][1, :]) > 1.0
    @test performance(r).feed[3] > 0
end

@testset "The default light product follows the key component" begin
    yl = M.light_product_composition([0.03, 0.15, 0.82]; component = 2)
    @test yl[2] == 1e-10 && sum(yl) ≈ 1 && yl[1] / yl[3] ≈ 0.03 / 0.82
    # CO2 is captured but is not the first component: the light product, and
    # the boundaries that take it, are stripped of CO2 rather than H2O
    gas = Gas(["H2O", "CO2", "N2"])
    function build(feed)
        local f = Flowsheet(gas, feed)
        add_bed!(f, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column(); ncells = 5)
        add_device!(f, :V_feed, :feed; at = (:Bed, :bottom))
        add_device!(f, :V_product, :product; at = (:Bed, :top))
        add_device!(f, :V_purge, :light_inlet; at = (:Bed, :top))
        add_device!(f, :V_vacuum, :extract; at = (:Bed, :bottom))
        return f
    end
    cycle = Cycle([Pressurise(5.0; source = :light_product), Adsorb(10.0; superficial_velocity = 0.37), Evacuate(10.0; to = 0.1bar)])
    fs = build(Feed(gas; composition = (H2O = 0.03, CO2 = 0.15, N2 = 0.82), pressure = 1bar))
    p = Process(fs, cycle; key_component = "CO2")
    yl = p.flowsheet.feed.light_product
    @test yl[2] == 1e-10 && yl[1] > 0.01 && sum(yl) ≈ 1
    @test M.boundary_condition(p.flowsheet, :V_purge, :inlet).composition == yl
    @test M.boundary_condition(p.flowsheet, :V_product, :outlet).composition == yl
    @test fs.feed.light_product[1] == 1e-10  # the flowsheet given is unchanged
    # A light product that is given is kept
    fs_given = build(Feed(gas; composition = (H2O = 0.03, CO2 = 0.15, N2 = 0.82), pressure = 1bar, light_product = (N2 = 1.0,)))
    @test Process(fs_given, cycle; key_component = "CO2").flowsheet.feed.light_product == fs_given.feed.light_product
end

@testset "Steps, roles and cycles" begin
    gas = M.haghpanah_2013_gas()
    feed = Feed(gas; composition = [0.15, 0.85], pressure = 1bar)
    a = load_adsorbent("zeolite_13x_haghpanah_2013")
    col = M.haghpanah_2013_column()
    fs = Flowsheet(gas, feed)
    add_bed!(fs, :Bed, a, col; ncells = 5)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_purge, :light_inlet; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    @test devices_with_role(fs, :extract) == [:V_vacuum]
    @test_throws MoccaInputError add_device!(fs, :V_x, :vacuum; at = (:Bed, :bottom))
    # A device that cannot be connected leaves the flowsheet as it was, so the
    # call can be corrected under the same name
    fs2 = deepcopy(fs)
    n_units, n_connections = length(M.unit_names(fs2)), length(M.connections(fs2))
    @test_throws MoccaInputError add_device!(fs2, :V_waste, :waste; at = (:Bed, :botom))
    @test_throws MoccaInputError add_device!(fs2, :V_waste, :link; from = (:Bed, :bottom), to = (:Bed, :middle))
    @test length(M.unit_names(fs2)) == n_units && length(M.connections(fs2)) == n_connections
    @test !haskey(fs2, :V_waste) && isnothing(M.device_role(fs2, :V_waste))
    @test !any(r -> r.name == :V_waste, keys(fs2.network.boundaries))
    add_device!(fs2, :V_waste, :waste; at = (:Bed, :bottom))
    @test devices_with_role(fs2, :waste) == [:V_waste]
    cycle = Cycle([Pressurise(10.0), Adsorb(20.0; superficial_velocity = 0.37), Blowdown(10.0; to = 0.2bar),
        Evacuate(20.0; to = 0.1bar), Purge(5.0; velocity = 0.1), Idle(5.0)])
    c = M.compile_cycle(cycle, fs)
    st = c.stages
    @test [s.name for s in st] == ["pressurisation", "adsorption", "blowdown", "evacuation", "purge", "idle"]
    # Pressures carry around the cycle: pressurisation starts where the purge
    # (at the evacuation pressure) leaves the bed
    @test st[1].settings[:V_feed].inlet.pressure.start == 0.1bar
    @test st[1].settings[:V_feed].inlet.pressure.stop == 1bar
    @test st[3].settings[:V_product].outlet.pressure.start == 1bar
    @test st[4].settings[:V_vacuum].outlet.pressure.start == 0.2bar
    @test st[2].settings[:V_feed].law.rate ≈ 0.37 * π * col.radius^2
    @test only(c.velocities[:Bed]) ≈ 1.0
    @test st[5].settings[:V_purge].law isa VolumetricFlow && st[5].settings[:V_vacuum].outlet.pressure == 0.1bar
    @test isempty(st[6].settings)
    @test st[3].settings[:V_product].law.conductance == port_conductance(fs[:Bed], 5)
    # Errors name the step
    err(c) = try M.compile_cycle(c, fs); nothing catch e; e end
    e = err(Cycle([Pressurise(10.0), Adsorb(10.0; velocity = 1.0), Blowdown(10.0; to = 2bar), Evacuate(10.0; to = 0.1bar)]))
    @test e isa MoccaInputError && occursin("step 3 (blowdown)", e.location)
    e = err(Cycle([Pressurise(10.0), Adsorb(10.0), Evacuate(10.0; to = 0.1bar)]))
    @test e isa MoccaInputError && occursin("velocity", e.message)
    e = err(Cycle([Pressurise(10.0), Adsorb(10.0; velocity = 1.0), Evacuate(10.0; to = 0.1bar, through = :middle)]))
    @test e isa MoccaInputError && e.field == "through"
    @test err(Cycle([Adsorb(-1.0; velocity = 1.0)])) isa MoccaInputError

    # An extension: a step type defined outside Mocca (above), registered for files
    cyc =Cycle([HoldAtFeed(duration = 5.0), Adsorb(10.0; velocity = 1.0), Blowdown(10.0; to = 0.5bar), Evacuate(10.0; to = 0.1bar)])
    p = Process(fs, cyc; timestepping = TimeStepping(max_dt = 2.0, first_dt = 0.5))
    @test M.stages(p)[1].name == "hold"
    @test M.stages(p)[1].settings[:V_feed].inlet.pressure == 1bar
    q = read_process(write_process(joinpath(mktempdir(), "hold.json"), p))
    @test M.stages(q)[1].name == "hold" && M.cycle_steps(q.cycle)[1][2] isa HoldAtFeed
    @test isempty(validate(p))

    # Equalisation through a device given by name: it must join the two beds
    fs2 = two_bed_vsa_process(ncells = 5).flowsheet
    eq(device) = M.compile_cycle(Cycle([Equalise(5.0; beds = (:A, :B), device = device)]), fs2)
    @test only(eq(:V_eq).stages).settings[:V_eq].law isa M.LinearValve
    @test_throws MoccaInputError eq(:V_vacuum_A)
    @test_throws MoccaInputError eq(:Nowhere)
end

@testset "Processes, files and validation" begin
    p = four_stage_vsa_process(ncells = 10, timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05))
    add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.5bar), description = "evacuation pressure")
    add_variable!(p, :rho, (:parameter, :Bed, :AdsorbentDensity); bounds = (1000.0, 1200.0))
    p[:rho] = 1150.0
    p[:P_L] = 0.12bar
    @test M.stages(p)[4].settings[:V_vacuum].outlet.pressure.stop == 0.12bar
    @test M.stages(p)[1].settings[:V_feed].inlet.pressure.start == 0.12bar
    # The extract line sits at the lowest pressure of the cycle
    @test M.boundary_condition(p.flowsheet, :V_vacuum, :outlet).pressure == 0.12bar
    @test_throws MoccaInputError add_variable!(p, :bad, ("evacuation", :speed))
    @test_throws MoccaInputError add_variable!(p, :bad, ("evacuationn", :to))
    path = write_process(joinpath(mktempdir(), "p.json"), p)
    q = read_process(path)
    @test q[:P_L] == 0.12bar && q[:rho] == 1150.0 && q.parameters[:Bed][:AdsorbentDensity] == p.parameters[:Bed][:AdsorbentDensity]
    fa, = M.cycle_schedule(p)
    fb, = M.cycle_schedule(q)
    @test force_difference(fa, p.model, fb, q.model) == 0
    @test library_processes() == ["haghpanah_2013_four_stage_vsa", "krishnamurthy_2014_dual_adsorbent_vsa", "krishnamurthy_2014_lpp_vsa"]
    @test read_process("haghpanah_2013_four_stage_vsa") isa Process
    # A hand-written file with units in strings
    txt = """
    {"gas": {"components": ["CO2", "N2"]},
     "feed": {"composition": {"CO2": 0.15, "N2": 0.85}, "pressure": "1 bar", "temperature": "25 °C"},
     "adsorbents": {"13X": "zeolite_13x_haghpanah_2013"},
     "columns": {"c": {"length": "1 m", "radius": "14.45 cm", "outer_radius": "16.2 cm", "porosity": 0.37}},
     "units": [{"name": "Bed", "adsorbent": "13X", "column": "c", "ncells": 10}],
     "devices": [{"name": "V_feed", "role": "feed", "at": ["Bed", "bottom"]},
                 {"name": "V_product", "role": "product", "at": ["Bed", "top"]},
                 {"name": "V_vacuum", "role": "extract", "at": ["Bed", "bottom"]}],
     "cycle": {"steps": [{"type": "pressurise", "duration": "15 s"}, {"type": "adsorb", "duration": 15, "velocity": "100 cm/s"},
                         {"type": "blowdown", "duration": 30, "to": "0.2 bar"}, {"type": "evacuate", "duration": 40, "to": "100 mbar"}]}}
    """
    file = joinpath(mktempdir(), "hand.json")
    write(file, txt)
    h = read_process(file)
    @test M.stages(h)[4].settings[:V_vacuum].outlet.pressure.stop ≈ 0.1bar
    @test h.flowsheet.feed.temperature ≈ 298.15
    for (from, to) in ("\"evacuate\"" => "\"evacuat\"", "0.2 bar" => "0.2 barr", "\"role\": \"extract\"" => "\"role\": \"vacuum\"",
            "\"ncells\": 10" => "\"ncells\": 10, \"cells\": 3")
        write(file, replace(txt, from => to))
        @test_throws MoccaInputError read_process(file)
    end
    # Dispersion set for a velocity far from the adsorption step's is reported
    gas = M.haghpanah_2013_gas()
    fs = Flowsheet(gas, Feed(gas; composition = [0.15, 0.85], pressure = 1bar))
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column(); ncells = 5, dispersion_velocity = 5.0)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    v = validate(Process(fs, Cycle([Pressurise(5.0; from = 0.5bar), Adsorb(5.0; velocity = 1.0)])))
    @test any(e -> e.field == "dispersion_velocity", v)
    @test any(e -> occursin("extract", e.message), v)
end

# Header and rows of a CSV file, with quoted fields
function read_csv(path; delim = ',')
    rows = map(readlines(path)) do line
        fields, buf, quoted = String[], IOBuffer(), false
        for c in line
            if c == '"'
                quoted = !quoted
            elseif c == delim && !quoted
                push!(fields, String(take!(buf)))
            else
                write(buf, c)
            end
        end
        push!(fields, String(take!(buf)))
    end
    return first(rows), rows[2:end]
end

@testset "Results written to CSV files" begin
    ncells = 10
    p = four_stage_vsa_process(ncells = ncells, timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05))
    r = simulate(p; cycles = 2)
    perf = performance(r)
    dir = mktempdir()
    paths = write_results(dir, r)
    @test sort(basename.(paths)) == sort(["summary.csv", "streams.csv", "stage_streams.csv", "stream_history.csv",
        "bed_history.csv", "bed_profiles.csv", "process.json"])
    @test read_process(joinpath(dir, "process.json")) isa Process
    col(header, name) = findfirst(==(name), header)
    num(x) = parse(Float64, x)

    h, rows = read_csv(joinpath(dir, "summary.csv"))
    @test h == ["Section", "Quantity", "Value", "Unit"]
    value(q, unit) = num(only(row[3] for row in rows if row[2] == q && row[4] == unit))
    @test value("CO2 purity in extract", "mol%") ≈ 100 * perf.purity rtol = 1e-6
    @test value("CO2 recovery", "%") ≈ 100 * perf.recovery rtol = 1e-6
    @test value("Specific energy", "kWh/t CO2") ≈ perf.specific_energy rtol = 1e-6
    @test value("Pressure", "bara") ≈ p.flowsheet.feed.pressure / 1e5
    # Away from the cyclic steady state, what the beds gain closes the balance
    @test abs(value("Total imbalance", "% of in")) < 1e-4
    @test abs(value("CO2 accumulated in beds", "mol/cycle")) > 1e-3

    h, streams = read_csv(joinpath(dir, "streams.csv"))
    @test [s[col(h, "Device")] for s in streams] == ["V_feed", "V_product", "V_vacuum"]
    extract = only(s for s in streams if s[col(h, "Role")] == "extract")
    @test extract[col(h, "From")] == "Bed bottom" && extract[col(h, "To")] == "boundary"
    @test num(extract[col(h, "CO2 [mol/cycle]")]) ≈ perf.extract[1] rtol = 1e-6
    @test num(extract[col(h, "CO2 [mol%]")]) ≈ 100 * perf.purity rtol = 1e-6

    # The stages add up to the cycle
    hs, stage_rows = read_csv(joinpath(dir, "stage_streams.csv"))
    @test unique(s[col(hs, "Stage")] for s in stage_rows) == [st.name for st in M.stages(p)]
    for s in streams
        d = s[col(h, "Device")]
        @test sum(num(t[col(hs, "N2 [mol/stage]")]) for t in stage_rows if t[col(hs, "Device")] == d) ≈
            num(s[col(h, "N2 [mol/cycle]")]) rtol = 1e-6
    end

    nsteps = length(r.timesteps)
    h, history = read_csv(joinpath(dir, "stream_history.csv"))
    @test length(history) == 3 * nsteps
    @test sort(unique(s[col(h, "Cycle")] for s in history)) == ["1", "2"]
    @test num(history[end][col(h, "Time [s]")]) ≈ 2 * M.cycle_time(p)
    # Nothing flows through a closed device, so it has no composition
    closed = [s for s in history if s[col(h, "Stage")] == "evacuation" && s[col(h, "Device")] == "V_feed"]
    @test all(s -> num(s[col(h, "Flow [mol/s]")]) == 0 && isempty(s[col(h, "CO2 [mol%]")]), closed)

    h, bed = read_csv(joinpath(dir, "bed_history.csv"))
    @test length(bed) == 2 * nsteps
    last_top = bed[end]
    @test last_top[col(h, "End")] == "top"
    @test num(last_top[col(h, "Pressure [bara]")]) ≈ r.states[end][:Bed][:Pressure][end] / 1e5 rtol = 1e-6
    @test num(last_top[col(h, "q CO2 [mol/kg]")]) ≈
        r.states[end][:Bed][:AdsorbedConcentration][1, end] / first(p.parameters[:Bed][:AdsorbentDensity]) rtol = 1e-6

    h, profiles = read_csv(joinpath(dir, "bed_profiles.csv"))
    @test length(profiles) == (length(M.stages(p)) + 1) * ncells
    @test profiles[1][col(h, "Snapshot")] == "start of cycle" && profiles[end][col(h, "Snapshot")] == "end of evacuation"
    @test num(profiles[end][col(h, "Temperature [degC]")]) ≈ r.states[end][:Bed][:Temperature][end] - 273.15 rtol = 1e-6

    # For spreadsheets with a decimal comma
    dir2 = mktempdir()
    write_results(dir2, r; delim = ';', decimal = ',')
    h2, streams2 = read_csv(joinpath(dir2, "streams.csv"); delim = ';')
    @test h2 == first(read_csv(joinpath(dir, "streams.csv")))
    @test num(replace(streams2[3][col(h2, "CO2 [mol/cycle]")], ',' => '.')) ≈ perf.extract[1] rtol = 1e-6
    @test_throws MoccaInputError write_results(dir2, r; delim = ',', decimal = ',')
    # Net amounts of components that went opposite ways have no composition
    @test all(isnan, M.Interface._fractions_percent([-0.5, 0.02]))
    @test M.Interface._fractions_percent([-1.0, -3.0]) ≈ [25.0, 75.0]
end

@testset "A flowsheet holds a unit network" begin
    gas = M.haghpanah_2013_gas()
    feed = Feed(gas; composition = [0.15, 0.85], pressure = 1bar)
    fs = Flowsheet(gas, feed)
    @test fs isa M.AbstractFlowsheet && M.network(fs) isa UnitNetwork
    @test fieldtype(Flowsheet, :gas) == Gas && fieldtype(Flowsheet, :feed) == Feed
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column(); ncells = 5)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    @test M.unit_names(fs) == M.unit_names(M.network(fs)) == [:Bed, :V_feed, :V_vacuum]
    @test M.holdup_units(fs) == M.beds(fs) == [:Bed]
    @test M.attachments(fs, :V_vacuum) == [M.StreamConnection(M.PortRef(:Bed, :bottom), M.PortRef(:V_vacuum, :inlet))]
    # The lower-level functions work on the flowsheet
    @test M.setup_flowsheet_model(fs) isa Jutul.MultiModel
    # Role boundaries are defaults until set
    r = M.PortRef(:V_vacuum, :outlet)
    @test r in fs.default_boundaries
    set_boundary!(fs, :V_vacuum, :outlet, PortCondition(pressure = 0.2bar, temperature = 298.15, composition = feed.composition))
    @test !(r in fs.default_boundaries) && M.boundary_condition(fs, :V_vacuum, :outlet).pressure == 0.2bar
    # Building on a bare network needs a Flowsheet's gas and feed
    @test_throws MethodError add_bed!(UnitNetwork(), :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column())
end

@testset "A process owns its flowsheet and cycle" begin
    gas = M.haghpanah_2013_gas()
    feed = Feed(gas; composition = [0.15, 0.85], pressure = 1bar)
    fs = Flowsheet(gas, feed)
    add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"), M.haghpanah_2013_column(); ncells = 5)
    add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs, :V_product, :product; at = (:Bed, :top))
    add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    cycle = Cycle([Pressurise(5.0), Adsorb(5.0; velocity = 0.5), Blowdown(5.0; to = 0.3bar), Evacuate(5.0; to = 0.1bar)])
    bed, boundary, spec = fs[:Bed], M.boundary_condition(fs, :V_vacuum, :outlet), deepcopy(fs.specs[:Bed])
    p = Process(fs, cycle)
    q = Process(fs, cycle)
    # Building them leaves the flowsheet and cycle as they were
    @test fs[:Bed] === bed && M.boundary_condition(fs, :V_vacuum, :outlet) == boundary
    @test fs.specs[:Bed].velocity_used == spec.velocity_used
    # The model of each process is built from its own units
    for r in (p, q)
        @test all(r.model.models[u] === r.flowsheet[u] for u in M.unit_names(r.flowsheet))
    end
    # Changing one process changes neither the other nor what they came from
    add_variable!(p, :P_L, ("evacuation", :to))
    add_variable!(p, :P_H, (:feed, :pressure))
    p[:P_L] = 0.12bar
    p[:P_H] = 1.2bar
    @test M.cycle_steps(cycle)[4][2].to == 0.1bar && feed.pressure == 1bar
    @test M.cycle_steps(q.cycle)[4][2].to == 0.1bar && q.flowsheet.feed.pressure == 1bar
    @test M.stages(q)[4].settings[:V_vacuum].outlet.pressure.stop == 0.1bar
    @test M.stages(p)[4].settings[:V_vacuum].outlet.pressure.stop == 0.12bar
    # A new feed temperature reaches the stages at once
    add_variable!(p, :T_feed, (:feed, :temperature))
    p[:T_feed] = 310.0
    @test M.stages(p)[1].settings[:V_feed].inlet.temperature == 310.0
    @test M.stages(p)[2].settings[:V_feed].inlet.temperature == 310.0
end

@testset "Axial dispersion is a parameter" begin
    gas = M.haghpanah_2013_gas()
    a = load_adsorbent("zeolite_13x_haghpanah_2013")
    function process(; v = 0.5, kwarg...)
        fs = Flowsheet(gas, Feed(gas; composition = [0.15, 0.85], pressure = 1bar))
        add_bed!(fs, :Bed, a, M.haghpanah_2013_column(); ncells = 5, kwarg...)
        add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
        add_device!(fs, :V_product, :product; at = (:Bed, :top))
        add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
        return Process(fs, Cycle([Pressurise(5.0), Adsorb(5.0; velocity = v), Evacuate(5.0; to = 0.1bar)]))
    end
    # Set from the adsorption velocity, without rebuilding the bed
    D(v) = M.compute_dispersion(M.axial_diffusivity(gas), v, a.particle_diameter)
    p = process(v = 0.5)
    @test only(p.parameters[:Bed][:AxialDispersion]) ≈ D(0.5)
    @test p.flowsheet.specs[:Bed].velocity_used == 0.5
    # A fixed velocity or dispersion is kept
    @test only(process(v = 0.5; dispersion_velocity = 2.0).parameters[:Bed][:AxialDispersion]) ≈ D(2.0)
    @test only(process(v = 0.5; dispersion = 1e-3).parameters[:Bed][:AxialDispersion]) == 1e-3
    # Same flux as the transmissibility of the old DiffusionTransmissibilities,
    # porosity·D through the faces
    bed = p.flowsheet[:Bed]
    g = Jutul.physical_representation(bed.data_domain)
    T_old = Jutul.compute_face_trans(g, bed.data_domain[:porosity] .* D(0.5))
    @test D(0.5) .* p.parameters[:Bed][:DispersionTransmissibilities] ≈ T_old
    # And a design variable can set it
    add_variable!(p, :D_ax, (:parameter, :Bed, :AxialDispersion))
    p[:D_ax] = 2e-3
    @test only(p.parameters[:Bed][:AxialDispersion]) == 2e-3
end

@testset "Design variable gradients and optimisation" begin
    p = four_stage_vsa_process(ncells = 8, durations = (10.0, 15.0, 20.0, 25.0), P_H = 1.1bar,
        timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05))
    add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.3bar))
    add_variable!(p, :t_ads, ("adsorption", :duration); bounds = (10.0, 30.0))
    add_variable!(p, :P_H, (:feed, :pressure); bounds = (1.05bar, 1.5bar))
    io = IOBuffer()
    g = performance_gradient(p; progress = io)
    @test g.steady_state.converged
    out = String(take!(io))
    @test startswith(out, "Cyclic steady state of four-stage VSA: 3 plain cycles, then Newton's method to 1.0e-09\n")
    @test contains(out, "  Cycles ") && contains(out, "  Newton ")
    @test contains(out, "done, on to Newton") && contains(out, "converged")
    # Every cycle simulated: one per line, and any step tried and rejected
    n_lines = count(contains(r"(cycle|iteration) +\d+  "), split(out, '\n'))
    @test parse(Int, match(r"Converged after (\d+) cycles simulated", out)[1]) >= n_lines
    value() = performance(cyclic_steady_state(p; tol = 1e-11))
    # Newton–Krylov from the last steady state, its Jacobian as preconditioner
    k = cyclic_steady_state(p; tol = 1e-11, krylov = true, progress = false)
    @test k.converged
    @test performance(k).recovery ≈ g.metrics.recovery rtol = 1e-8
    for n in (:P_L, :t_ads, :P_H)
        x0 = p[n]
        h = 1e-4 * x0
        p[n] = x0 + h; up = value()
        p[n] = x0 - h; down = value()
        p[n] = x0
        for k in (:purity, :recovery, :specific_energy, :productivity)
            @test g.gradients[k][n] ≈ (up[k] - down[k]) / 2h rtol = 1e-4
        end
    end
    r = optimise!(p; variables = [:P_L, :t_ads], maximise = :recovery, at_least = (purity = 0.5,), max_it = 2)
    @test r.evaluations >= 2
    @test 0.05bar <= r.values.P_L <= 0.3bar && 10.0 <= r.values.t_ads <= 30.0
    @test p[:P_L] == r.values.P_L
    @test r.metrics.recovery >= first(r.history).metrics.recovery - 1e-9
end

@testset "A fine grid starts from the steady state of a coarser one" begin
    I = M.Interface
    p = four_stage_vsa_process(ncells = 20)
    # Half the cells and the same steps; the process itself is left as it was
    q = I._coarse_copy(p, 5)
    @test q.flowsheet.specs[:Bed].ncells == 10 && p.flowsheet.specs[:Bed].ncells == 20
    @test length(M.initial_state(q)[:Bed][:Pressure]) == 10
    @test I.cycle_schedule(q)[2] == I.cycle_schedule(p)[2]
    @test isnothing(I._coarse_copy(p, 20)) && isnothing(I._coarse_copy(p, nothing))
    @test [g.flowsheet.specs[:Bed].ncells for g in I._coarse_grids(p, 5)] == [5, 10]
    # Linear interpolation, cell centre to cell centre: a linear profile is
    # kept inside the bed, the end cells take the nearest coarse cell, and the
    # mole fractions still add up to one
    x = M.initial_state(q)
    x[:Bed][:Pressure] .= 1e5 .* (1 .+ ((1:10) .- 0.5) ./ 10)
    x[:Bed][:y][1, :] .= ((1:10) .- 0.5) ./ 10
    x[:Bed][:y][2, :] .= 1 .- x[:Bed][:y][1, :]
    xf = I._interpolate_state(p, x)
    @test xf[:Bed][:Pressure][2:19] ≈ 1e5 .* (1 .+ ((2:19) .- 0.5) ./ 20)
    @test xf[:Bed][:Pressure][[1, 20]] ≈ x[:Bed][:Pressure][[1, 10]]
    @test vec(sum(xf[:Bed][:y]; dims = 1)) ≈ ones(20)
    @test I._interpolate_state(q, x)[:Bed][:Pressure] == x[:Bed][:Pressure]

    # The coarsest grid from the initial state, each finer one by Newton's
    # method from the one before, to the same steady state as from the
    # initial state on the fine grid. Each grid decides how to solve Newton's
    # steps: here GMRES above 100 variables, so only on 20 cells (120)
    io = IOBuffer()
    r = cyclic_steady_state(p; coarse_cells = 5, krylov = 100, progress = io)
    @test r.converged
    out = String(take!(io))
    @test startswith(out, "Cyclic steady state of four-stage VSA: 3 plain cycles, then Newton's method on 5 cells to 1.0e-06; " *
        "from there Newton's method on 10 cells to 1.0e-06 and on 20 cells to 1.0e-09\n")
    @test contains(out, "10 cells: Newton's method from the steady state on 5 cells")
    @test contains(out, "20 cells: Newton's method from the steady state on 10 cells")
    grids = split(out, r"\n  \d+ cells: ")
    @test length(grids) == 4
    @test contains(grids[2], "Jacobian") && !contains(grids[2], "GMRES")
    @test contains(grids[3], "Jacobian") && !contains(grids[3], "GMRES")
    @test contains(grids[4], "GMRES") && !contains(grids[4], "Jacobian")
    @test_throws MoccaInputError cyclic_steady_state(p; krylov = :yes)
    cold = cyclic_steady_state(p; coarse_cells = nothing, warm_start = false, progress = false)
    @test cold.converged
    @test cycle_change(p.model, r.state0, cold.state0) < 1e-7
    @test performance(r).recovery ≈ performance(cold).recovery rtol = 1e-7
    # From its last steady state the process starts on its own grid
    io = IOBuffer()
    cyclic_steady_state(p; coarse_cells = 5, progress = io)
    @test startswith(String(take!(io)), "Cyclic steady state of four-stage VSA: Newton's method to 1.0e-09 from the last steady state\n")
end
