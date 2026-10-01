# `time_scale` sets the TimeScale parameter for the step; it is only used by
# the adjoints for stage durations
function Jutul.setup_forces(model::AdsorptionModel; bc = nothing, time_scale = nothing)
    return (bc = bc, time_scale = time_scale)
end

"""
    compute_column_face_area(model, state)

Flow cross-section of the bed [m²]. Read from the `BedCrossSectionArea`
parameter when `state` carries parameters, otherwise from the column radius in
the data domain.
"""
function compute_column_face_area(model::AdsorptionModel, state)
    if _state_has(state, :BedCrossSectionArea)
        return state.BedCrossSectionArea[1]
    else
        r_in = first(model.data_domain[:r_in, Unit()])
        return π * r_in^2
    end
end

"""
    calc_bc_trans(model, state, cell)

Transmissibility [m³] of the half-cell connection between `cell` and the bed
end it touches.
"""
function calc_bc_trans(model::AdsorptionModel, state, cell)
    k = _state_or_domain(state, model, :Permeability, :permeability)[cell]
    dx = _state_or_domain(state, model, :CellDx, :dx)[cell] / 2.0
    A = compute_column_face_area(model, state)
    return k * A / dx
end

include("bc_pressurisation.jl")
include("bc_adsorption.jl")
include("bc_blowdown.jl")
include("bc_evacuation.jl")