struct Column <: Jutul.JutulEntity end

# Overload JutulDarcy functions

JutulDarcy.component_names(sys::AdsorptionSystem) = sys.component_names
JutulDarcy.number_of_components(sys::AdsorptionSystem) = length(sys.component_names)
JutulDarcy.has_other_phase(::AdsorptionSystem) = false
JutulDarcy.phase_names(::AdsorptionSystem) = ["gas"]
JutulDarcy.number_of_phases(::AdsorptionSystem) = 1

# `AdsorptionSystem` is not a `JutulDarcy.MultiPhaseSystem` in the method table, so without this
# forward the default `Jutul.discretize_domain(::DataDomain, ...)` would wrap the mesh only.
# Delegate to the `MultiPhaseSystem` implementation (TPFA + entity propagation from `DataDomain`).
function Jutul.discretize_domain(d::Jutul.DataDomain, system::AdsorptionSystem, ::Val{:default}; kwarg...)
    return invoke(
        Jutul.discretize_domain,
        Tuple{Jutul.DataDomain, JutulDarcy.MultiPhaseSystem, Val{:default}},
        d, system, Val(:default);
        kwarg...,
    )
end

# Specific functions needed for Mocca

"""
    column_mesh(ncells, L, r_in)

Create a 1D `CartesianMesh` for a cylindrical adsorption column.
The cross-section is a square with the same area as the circular column
(side length `√(π r_in²)`).
"""
function column_mesh(ncells, L, r_in)
    dx = sqrt(π * r_in^2)
    return Jutul.CartesianMesh((ncells, 1, 1), (L, dx, dx))
end

function compute_permeability(porosity, d_p)
    return 4 / 150 * ((porosity / (1 - porosity))^2) * (d_p / 2)^2 * porosity
end

function compute_dispersion(D_m, V0_inter, d_p)
    return 0.7 * D_m + 0.5 * V0_inter * d_p
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
