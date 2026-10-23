# Writes the steady cycle of the four-step VSA of Haghpanah et al. (2013) to
# docs/src/assets/home_cycle.json, for the animation on the home page. The
# file is committed, so the docs build does not run this. Run it from the
# repository root after a change to the model or the library process:
#
#     julia --project=. docs/home_cycle_data.jl

using Mocca
using Mocca.Units
using Printf

process = Mocca.read_process("haghpanah_2013_four_stage_vsa")
process[:P_L] = 0.05bar
result = Mocca.cyclic_steady_state(process)
result.converged || error("the cyclic steady state did not converge")
perf = Mocca.performance(result)

bed = process.flowsheet[:Bed]
z = Mocca.axial_positions(bed)
z = z ./ maximum(z)
t = cumsum(result.timesteps)
t_cycle = t[end]

# Frames at even times through the cycle, each the state at the end of the
# first time step that reaches it; the page interpolates between them
nframes = 160
frames = unique([findfirst(>=(t_cycle * (k - 0.5) / nframes), t) for k in 1:nframes])
bed_states = [result.states[k][:Bed] for k in frames]
q = [Mocca.axial_profile(bed, s[:AdsorbedConcentration])[1, :] for s in bed_states]
q_max = maximum(maximum.(q))

num(x) = @sprintf("%.4g", x)
vec_json(v) = "[" * join(num.(v), ",") * "]"

stages = Mocca.stages(process)

# The pressure set at the valve that is open at the bed's end in each stage:
# the feed inlet while the bed is pressurised and fed, and the outlet while it
# is blown down and evacuated. It follows the ramps of the steps, without the
# pressure drop that the feed flow adds inside the bed.
stage_ends = cumsum([st.duration for st in stages])
open_valve = Dict("pressurisation" => (:V_feed, :InletPressure), "adsorption" => (:V_feed, :InletPressure),
    "blowdown" => (:V_product, :OutletPressure), "evacuation" => (:V_vacuum, :OutletPressure))
function valve_pressure(k)
    st = stages[something(findfirst(>=(t[k] - 1e-9), stage_ends), length(stages))]
    device, port = open_valve[st.name]
    return result.states[k][device][port][1]
end

open(joinpath(@__DIR__, "src", "assets", "home_cycle.json"), "w") do io
    println(io, "{")
    println(io, "  \"process\": \"haghpanah_2013_four_stage_vsa\",")
    println(io, "  \"P_L_bar\": ", num(process[:P_L] / bar), ",")
    println(io, "  \"purity\": ", num(perf.purity), ",")
    println(io, "  \"recovery\": ", num(perf.recovery), ",")
    println(io, "  \"specific_energy_kwh_per_t\": ", num(perf.specific_energy), ",")
    println(io, "  \"cycle_time_s\": ", num(t_cycle), ",")
    println(io, "  \"stages\": [", join(["{\"name\": \"$(st.name)\", \"duration\": $(num(st.duration))}" for st in stages], ", "), "],")
    println(io, "  \"z\": ", vec_json(z), ",")
    println(io, "  \"t\": ", vec_json(t[frames]), ",")
    println(io, "  \"pressure_valve_bar\": ", vec_json([valve_pressure(k) / bar for k in frames]), ",")
    println(io, "  \"loading\": [", join([vec_json(qk ./ q_max) for qk in q], ",\n    "), "],")
    println(io, "  \"y_co2\": [", join([vec_json(Mocca.axial_profile(bed, s[:y])[1, :]) for s in bed_states], ",\n    "), "]")
    println(io, "}")
end
println("purity $(perf.purity), recovery $(perf.recovery), cycle $(t_cycle) s, $(length(stages)) stages")
