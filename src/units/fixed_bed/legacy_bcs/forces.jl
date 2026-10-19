# `time_scale` sets the TimeScale parameter for the step; it is only used by
# the adjoints for stage durations
function Jutul.setup_forces(model::AdsorptionModel; bc = nothing, time_scale = nothing)
    return (bc = bc, time_scale = time_scale)
end

include("bc_pressurisation.jl")
include("bc_adsorption.jl")
include("bc_blowdown.jl")
include("bc_evacuation.jl")