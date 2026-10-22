# Beds meshed in the (z, r) plane of the column: rings of the cylinder, with
# dispersion and conduction across the radius, and the wall around the outer
# ring.

using Mocca.Units

@testset "Axisymmetric column geometry" begin
    c = haghpanah_constants()
    radii = [0.0, 0.05, 0.11, c.r_in]
    bed = haghpanah_radial_bed(c; ncells = 5, radial_cells = 3, radii = radii)
    d = bed.data_domain
    L = sum(d[:dx][1:5])
    ring = [π * (radii[j + 1]^2 - radii[j]^2) for j in 1:3]
    @test sum(d[:volumes]) ≈ π * c.r_in^2 * L
    @test d[:volumes] ≈ repeat(ring, inner = 5) .* (L / 5)
    # A face across the axis is its ring, and one at a radius its cylinder
    N, a, ax = d[:neighbors], d[:areas], Mocca.Balances._axial_faces(d)
    for f in axes(N, 2)
        j = cld(N[1, f], 5)
        @test a[f] ≈ (ax[f] ? ring[j] : 2π * radii[j + 1] * L / 5)
    end
    # The ports are the end rings, sharing the stream as their areas do, and
    # conduct as the ends of the column
    @test Mocca.ports(bed) == (bottom = [1, 6, 11], top = [5, 10, 15])
    @test Mocca.port_weights(bed, [1, 6, 11]) ≈ ring ./ sum(ring)
    @test Mocca.port_conductance(bed, :top) ≈ Mocca.port_conductance(haghpanah_bed(c; ncells = 5), :top)
    # The wall surrounds the outer ring only, and conducts along it
    prm = Mocca.setup_process_parameters(bed)
    @test Mocca.wall_cells(bed) == 11:15
    @test prm[:WallCells] == [zeros(10); ones(5)]
    @test sum(prm[:WallAreaIn]) ≈ 2π * c.r_in * L && all(iszero, prm[:WallAreaIn][1:10])
    @test findall(>(0), prm[:WallTransmissibilities]) == [f for f in axes(N, 2) if ax[f] && N[1, f] > 10]
    # Dispersion along the axis and across it are separate parameters
    @test all(iszero, prm[:DispersionTransmissibilities][.!ax]) && all(iszero, prm[:RadialDispersionTransmissibilities][ax])
    @test prm[:RadialDispersion] == prm[:AxialDispersion]
    # One ring is the 1D column
    @test Mocca.column_mesh(5, L, c.r_in; radial_cells = 1).dims == (5, 1, 1)
    @test_throws ArgumentError Mocca.column_mesh(5, L, c.r_in; radial_cells = 3, radii = [0.0, 0.1, c.r_in])
    @test !haskey(Mocca.setup_process_parameters(haghpanah_bed(c; ncells = 5)), :RadialDispersion)
end

@testset "Conduction across the radius" begin
    # On uniform rings, the two-point flux of T = r² is exact: its divergence
    # is 4λ over each cell inside the outer ring, which has no face outside
    c = haghpanah_constants()
    nz, nr, λ_r = 3, 5, 0.7
    bed = haghpanah_radial_bed(c; ncells = nz, radial_cells = nr, λ_r = λ_r, thermal = Mocca.Adiabatic())
    d = bed.data_domain
    prm = Mocca.setup_process_parameters(bed)
    N, ax = d[:neighbors], Mocca.Balances._axial_faces(d)
    T = d[:cell_centroids][2, :] .^ 2
    for (K, k) in ((prm[:ThermalConductivities], λ_r), (prm[:RadialDispersionTransmissibilities], c.Φ))
        div = zeros(nz * nr)
        for f in axes(N, 2)
            ax[f] && continue
            l, r = N[:, f]
            div[l] += K[f] * (T[r] - T[l])
            div[r] += K[f] * (T[l] - T[r])
        end
        inside = 1:(nz * (nr - 1))
        @test div[inside] ≈ 4k .* d[:volumes][inside]
    end
end

# The bed and the column are compared on the same time steps: a short first
# step keeps Newton's method from cutting the steps at the start of
# pressurisation, where round-off decides between them
const _GENTLE_START = 0.1

@testset "A bed without radial gradients is the column" begin
    # Adiabatic and fed evenly across its end, nothing varies across the radius
    one = haghpanah_vsa(; ncells = 10, thermal = Mocca.Adiabatic(), first_dt = _GENTLE_START)
    r1 = simulate(one)
    for radii in (nothing, [0.0, 0.03, 0.07, 0.12, haghpanah_constants().r_in])
        p = haghpanah_vsa(; ncells = 10, thermal = Mocca.Adiabatic(), first_dt = _GENTLE_START,
            bed_model = haghpanah_radial_bed(; ncells = 10, radial_cells = 4, radii = radii, thermal = Mocca.Adiabatic()))
        r = simulate(p)
        @test r.timesteps ≈ r1.timesteps
        bed = p.flowsheet[:Bed]
        s, s1 = r.states[end][:Bed], r1.states[end][:Bed]
        for k in (:Pressure, :Temperature, :y, :AdsorbedConcentration)
            @test Mocca.axial_profile(bed, s[k]) ≈ s1[k] rtol = 1e-6
        end
        @test s[:y][:, 1:10] ≈ s[:y][:, 31:40] rtol = 1e-8
        for dev in (:V_feed, :V_product, :V_vacuum)
            @test Mocca.stream_totals(r.states, r.timesteps, dev) ≈ Mocca.stream_totals(r1.states, r1.timesteps, dev) rtol = 1e-6
        end
    end
end

@testset "Radial gradients and the wall" begin
    one = haghpanah_vsa(; ncells = 10, first_dt = _GENTLE_START)
    r1 = simulate(one)
    function run(λ_r)
        p = haghpanah_vsa(; ncells = 10, first_dt = _GENTLE_START,
            bed_model = haghpanah_radial_bed(; ncells = 10, radial_cells = 4, λ_r = λ_r))
        return p, simulate(p)
    end
    # Conduction across the radius much faster than to the wall makes the
    # column, whose wall takes heat from the whole cross-section
    p, r = run(1e3)
    @test r.timesteps ≈ r1.timesteps
    bed = p.flowsheet[:Bed]
    s, s1 = r.states[end][:Bed], r1.states[end][:Bed]
    @test maximum(abs.(Mocca.axial_profile(bed, s[:Temperature]) .- s1[:Temperature])) < 1e-3
    @test s[:WallTemperature][31:40] ≈ s1[:WallTemperature] atol = 1e-3
    @test performance(r).recovery ≈ performance(r1).recovery rtol = 1e-5
    # Slower, the heat of adsorption leaves the centre more slowly than the
    # outer ring, which loses it to the wall
    p, r = run(0.5)
    ΔT = maximum(maximum(st[:Bed][:Temperature][1:10] .- st[:Bed][:Temperature][31:40]) for st in r.states)
    @test ΔT > 0.3
    # Inside the outer ring the wall temperature follows the gas
    s = r.states[end][:Bed]
    @test s[:WallTemperature][1:30] ≈ s[:Temperature][1:30]
    # The bed ends in the results report the wall of the outer ring
    rows = Mocca.Interface._bed_history_table(r).rows
    top = rows[findlast(row -> row[6] == "top", rows)]
    @test top[9] ≈ s[:WallTemperature][40] - 273.15
end

@testset "Rings that differ conserve moles" begin
    p = haghpanah_vsa(; ncells = 10, bed_model = haghpanah_radial_bed(; ncells = 10, radial_cells = 3, porosity = [0.3, 0.37, 0.45]))
    state0, prm = Mocca.initial_state(p), p.parameters
    r = simulate(p; state0 = state0)
    bed = p.flowsheet[:Bed]
    @test maximum(abs(s[:Bed][:Pressure][1] - s[:Bed][:Pressure][21]) for s in r.states) > 1.0
    held(s) = Mocca.component_inventory(bed, s[:Bed], prm[:Bed])
    n_in = Mocca.stream_totals(r.states, r.timesteps, :V_feed)
    n_out = Mocca.stream_totals(r.states, r.timesteps, :V_product) .+ Mocca.stream_totals(r.states, r.timesteps, :V_vacuum)
    @test held(r.states[end]) .- held(state0) ≈ n_in .- n_out rtol = 1e-7 atol = 1e-8 * sum(n_in)
end

# A column with the heat transfer a bed meshed across its radius needs; the
# values are for testing, not from a source
_radial_column(; kw...) = Column(; length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    axial_conductivity = 0.0903, radial_conductivity = 0.5, wall_heat_transfer = 30.0, kw...)

@testset "Beds meshed across their radius in a flowsheet" begin
    gas = Mocca.haghpanah_2013_gas()
    a = load_adsorbent("zeolite_13x_haghpanah_2013")
    fs = Flowsheet(gas, Feed(gas; composition = [0.15, 0.85], pressure = 1bar))
    # The column must give the heat transfer across the radius
    @test_throws MoccaInputError add_bed!(fs, :B, a, Mocca.haghpanah_2013_column(); ncells = 4, radial_cells = 3)
    @test_throws MoccaInputError add_bed!(fs, :B, a, _radial_column(wall_heat_transfer = nothing); ncells = 4, radial_cells = 3)
    add_bed!(fs, :B, a, _radial_column(wall_heat_transfer = nothing); ncells = 4, radial_cells = 3, thermal = Adiabatic())
    @test_throws MoccaInputError add_bed!(fs, :C, a, _radial_column(); ncells = 4, radial_cells = 3, radii = [0.0, 0.1445])
    @test_throws MoccaInputError add_bed!(fs, :C, a, _radial_column(); ncells = 4, radii = [0.0, 0.1445])
    @test_throws MoccaInputError _radial_column(radial_conductivity = -1.0)

    p = four_stage_vsa_process(ncells = 6, radial_cells = 3, column = _radial_column(), durations = (5.0, 5.0, 5.0, 5.0))
    bed = p.flowsheet[:Bed]
    @test length(Mocca.port_cells(bed, :bottom)) == 3
    prm = p.parameters[:Bed]
    # The wall coefficient is that of the gas next to the wall
    @test only(prm[:InnerHeatTransferCoeff]) == 30.0
    # The dispersion across the radius follows the adsorption velocity as the
    # dispersion along it does
    v, d_p = p.flowsheet.specs[:Bed].velocity_used, a.particle_diameter
    @test only(prm[:RadialDispersion]) ≈ Mocca.compute_radial_dispersion(Mocca.axial_diffusivity(gas), v, d_p)
    @test only(prm[:AxialDispersion]) ≈ Mocca.compute_dispersion(Mocca.axial_diffusivity(gas), v, d_p)
    @test only(prm[:RadialDispersion]) < only(prm[:AxialDispersion])
    @test contains(sprint(Mocca.print_flowsheet_diagram, p.flowsheet), "6 × 3 cells")
    # A design variable can set it
    add_variable!(p, :D_r, (:parameter, :Bed, :RadialDispersion))
    p[:D_r] = 3e-4
    @test only(p.parameters[:Bed][:RadialDispersion]) == 3e-4

    # Written to a file and read back, the bed is the same
    q = four_stage_vsa_process(ncells = 6, radial_cells = 3, column = _radial_column(), durations = (5.0, 5.0, 5.0, 5.0))
    d = Mocca.process_to_dict(q)
    @test d["units"][1]["radial_cells"] == 3
    @test d["columns"][d["units"][1]["column"]]["wall_heat_transfer"] == 30.0
    path = tempname() * ".json"
    write_process(path, q)
    q2 = read_process(path)
    s1, s2 = q.flowsheet.specs[:Bed], q2.flowsheet.specs[:Bed]
    @test s2.radial_cells == 3 && s2.column.radial_conductivity == 0.5 && s2.column.wall_heat_transfer == 30.0
    @test all(q2.parameters[:Bed][k] ≈ q.parameters[:Bed][k] for k in keys(q.parameters[:Bed]))
    @test performance(simulate(q2)).recovery ≈ performance(simulate(q)).recovery
    # Uneven rings and a fixed radial dispersion are kept too
    fs2 = Flowsheet(gas, Feed(gas; composition = [0.15, 0.85], pressure = 1bar))
    add_bed!(fs2, :Bed, a, _radial_column(); ncells = 4, radial_cells = 2, radii = [0.0, 0.1, 0.1445], radial_dispersion = 1e-4)
    add_device!(fs2, :V_feed, :feed; at = (:Bed, :bottom))
    add_device!(fs2, :V_product, :product; at = (:Bed, :top))
    p2 = Process(fs2, Cycle([Pressurise(5.0), Adsorb(5.0; velocity = 0.5)]))
    @test only(p2.parameters[:Bed][:RadialDispersion]) == 1e-4
    write_process(path, p2)
    s = read_process(path).flowsheet.specs[:Bed]
    @test s.radii == [0.0, 0.1, 0.1445] && s.radial_dispersion == 1e-4
    # A bed of one radial cell writes nothing new
    @test !haskey(Mocca.process_to_dict(four_stage_vsa_process(ncells = 6))["units"][1], "radial_cells")
end

@testset "Gradient with respect to the radial dispersion" begin
    p = four_stage_vsa_process(ncells = 4, radial_cells = 2, column = _radial_column(), durations = (10.0, 15.0, 20.0, 25.0),
        P_H = 1.1bar, timestepping = TimeStepping(max_dt = 1.0, first_dt = 0.05))
    add_variable!(p, :D_r, (:parameter, :Bed, :RadialDispersion))
    g = performance_gradient(p; progress = false)
    @test g.steady_state.converged
    # The radial dispersion changes the performance little, so the finite
    # difference takes a large step to rise above the tolerance of the steady
    # state
    value() = performance(cyclic_steady_state(p; tol = 1e-11, progress = false))
    x0 = p[:D_r]
    h = 1e-2 * x0
    p[:D_r] = x0 + h; up = value()
    p[:D_r] = x0 - h; down = value()
    p[:D_r] = x0
    for k in (:purity, :recovery, :specific_energy, :productivity)
        @test g.gradients[k][:D_r] ≈ (up[k] - down[k]) / 2h rtol = 1e-3
    end
end

@testset "The steady state of a radial bed starts from one radial cell" begin
    I = Mocca.Interface
    c = _radial_column()
    p = four_stage_vsa_process(ncells = 20, radial_cells = 3, column = c)
    add_variable!(p, :D_r, (:parameter, :Bed, :RadialDispersion))
    p[:D_r] = 3e-4
    # One radial cell on as many cells along the column, then halving as for
    # any bed; the settings of parameters only the radial bed has are left out
    grids = I._coarse_grids(p, 5)
    @test [(g.flowsheet.specs[:Bed].ncells, g.flowsheet.specs[:Bed].radial_cells) for g in grids] == [(5, 1), (10, 1), (20, 1)]
    q = last(grids)
    @test !haskey(q.parameters[:Bed], :RadialDispersion) && p.flowsheet.specs[:Bed].radial_cells == 3
    # The 1D bed has the wall coefficient that gives the radial bed's mean
    # temperature
    @test only(q.parameters[:Bed][:InnerHeatTransferCoeff]) ≈ 1 / (1 / c.wall_heat_transfer + c.radius / (4 * c.radial_conductivity))

    # Spread over the rings, with the wall temperature of the gas inside the
    # outer ring
    x = Mocca.initial_state(q)
    x[:Bed][:Pressure] .= 1e5 .* (1 .+ ((1:20) .- 0.5) ./ 20)
    x[:Bed][:Temperature] .= 300 .+ (1:20)
    x[:Bed][:WallTemperature] .= 290.0
    xf = I._interpolate_state(p, x, q)
    for j in 1:3
        @test xf[:Bed][:Pressure][(1:20) .+ 20(j - 1)] ≈ x[:Bed][:Pressure]
    end
    @test xf[:Bed][:WallTemperature][41:60] ≈ fill(290.0, 20)
    @test xf[:Bed][:WallTemperature][1:40] ≈ xf[:Bed][:Temperature][1:40]
    @test vec(sum(xf[:Bed][:y]; dims = 1)) ≈ ones(60)
    # From rings to more rings, a field linear along the column and across
    # the radius is kept inside the bed
    p6 = four_stage_vsa_process(ncells = 20, radial_cells = 6, column = c)
    p3 = four_stage_vsa_process(ncells = 10, radial_cells = 3, column = c)
    f(m) = [m.data_domain[:cell_centroids][1, k] + 3 * m.data_domain[:cell_centroids][2, k] for k in 1:Jutul.number_of_cells(m.domain)]
    x3 = Mocca.initial_state(p3)
    x3[:Bed][:Temperature] .= 300 .+ f(p3.flowsheet[:Bed])
    T6 = I._interpolate_state(p6, x3, p3)[:Bed][:Temperature]
    inside = [i + 20(j - 1) for j in 2:5 for i in 2:19]
    @test T6[inside] .- 300 ≈ f(p6.flowsheet[:Bed])[inside]

    # Through the 1D bed to the same steady state as from the initial state
    small() = four_stage_vsa_process(ncells = 8, radial_cells = 2, column = c)
    io = IOBuffer()
    r = cyclic_steady_state(small(); coarse_cells = 4, progress = io)
    out = String(take!(io))
    @test r.converged
    @test contains(out, "8 × 2 cells: Newton's method from the steady state on 8 cells")
    cold = cyclic_steady_state(small(); coarse_cells = nothing, progress = false)
    @test cold.converged
    @test Mocca.cycle_change(small().model, r.state0, cold.state0) < 1e-8
end
