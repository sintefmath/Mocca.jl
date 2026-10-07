# Plot the fronts written by purity_recovery.jl and energy_productivity.jl, as
# in Figures 10a and 12 of Haghpanah et al. (2013). Every evaluated design is
# shown faintly behind the front.
using CairoMakie

const RESULTS = joinpath(@__DIR__, "..", "..", "results", "haghpanah_2013_moo")

function read_csv(path)
    lines = readlines(path)
    header = Symbol.(split(lines[1], ","))
    rows = [parse.(Float64, split(l, ",")) for l in lines[2:end] if !isempty(l)]
    return Dict(h => [r[i] for r in rows] for (i, h) in enumerate(header))
end

files(pattern) = sort(filter(f -> occursin(pattern, f), readdir(RESULTS)))
label_of(f, key) = (m = match(Regex("$(key)_([0-9.]+)"), f); isnothing(m) ? f : m[1])

fig = Figure(size = (1100, 460))
ax1 = Axis(fig[1, 1], xlabel = "CO₂ purity [%]", ylabel = "CO₂ recovery [%]",
    title = "Purity–recovery (cf. Figure 10)")
for f in files(r"^p?vsa_PLmin_.*_front\.csv$")
    ev = joinpath(RESULTS, replace(f, "_front" => "_evaluations"))
    if isfile(ev)
        e = read_csv(ev)
        scatter!(ax1, 100e[:purity], 100e[:recovery], color = (:gray, 0.15), markersize = 4)
    end
    d = read_csv(joinpath(RESULTS, f))
    o = sortperm(d[:purity])
    kind = startswith(f, "pvsa") ? "PVSA" : "VSA"
    scatterlines!(ax1, 100d[:purity][o], 100d[:recovery][o], label = "$kind, P_L ≥ $(label_of(f, "PLmin")) bar")
end
hlines!(ax1, [90.0], color = :black, linestyle = :dash, linewidth = 0.8)
vlines!(ax1, [90.0], color = :black, linestyle = :dash, linewidth = 0.8)
axislegend(ax1, position = :lb)

ax2 = Axis(fig[1, 2], xlabel = "productivity [mol CO₂ m⁻³ adsorbent s⁻¹]", ylabel = "energy [kWh t⁻¹ CO₂]",
    title = "Energy–productivity at purity, recovery ≥ 90% (cf. Figure 12)")
for f in files(r"^p?vsa_PL_.*_front\.csv$")
    d = read_csv(joinpath(RESULTS, f))
    ok = (d[:purity] .>= 0.895) .& (d[:recovery] .>= 0.895)
    any(ok) || continue
    o = sortperm(d[:productivity][ok])
    kind = startswith(f, "pvsa") ? "PVSA" : "VSA"
    scatterlines!(ax2, d[:productivity][ok][o], d[:energy][ok][o], label = "$kind, P_L = $(label_of(f, "PL")) bar")
end
axislegend(ax2, position = :lt)
save(joinpath(RESULTS, "fronts.png"), fig)
println("saved ", joinpath(RESULTS, "fronts.png"))
