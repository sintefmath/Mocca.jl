# Beds meshed in the (z, r) plane of the column: rings of the cylinder, with
# dispersion and conduction across the radius, and the wall around the outer
# ring.

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
