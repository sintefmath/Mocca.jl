function Jutul.setup_forces(model::AdsorptionModel; bc = nothing)
    return (bc = bc,)
end

function compute_column_face_area(model::AdsorptionModel, state)
    g = Jutul.physical_representation(model.data_domain)::Jutul.CartesianMesh
    return (g.deltas[2] * g.deltas[3])
end

function calc_bc_trans(model::AdsorptionModel, state, cell)
    k = model.data_domain[:permeability][cell]
    r_in = first(model.data_domain[:r_in, Column()])
    dx = state.CellDx[cell] / 2.0
    A = (π * r_in^2)
    return k * A / dx
end

function calc_bc_wall_trans(model::AdsorptionModel, state, cell)
    k = state.WallConductivity[1]
    dx = state.CellDx[cell] / 2.0
    A = state.WallCrossSectionArea[1]
    return k * A / dx
end

include("bc_pressurisation.jl")
include("bc_adsorption.jl")
include("bc_blowdown.jl")
include("bc_evacuation.jl")