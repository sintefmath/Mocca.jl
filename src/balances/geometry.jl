# ------------------------------------------------------------------------------
# Axisymmetric columns
# ------------------------------------------------------------------------------

# Volumes and face areas of the cells of the (z, r) plane of a column (see
# column_mesh) as rings of the cylinder: a cell from r_i to r_o and Δz long
# holds π(r_o² − r_i²)·Δz, a face across the axis is the ring's cross-section,
# and a face at radius r is 2π·r·Δz. The centroids stay in the plane, so the
# two-point fluxes run over distances along z and r.
function _axisymmetric_geometry!(domain, mesh, r_in)
    nz, nr, nk = Jutul.grid_dims_ijk(mesh)
    nk == 1 || throw(ArgumentError("an axisymmetric column has a mesh of one cell in z, got $nk"))
    radii = [0.0; cumsum(collect(mesh.deltas[2]))]
    last(radii) ≈ r_in || throw(ArgumentError("the mesh reaches r = $(last(radii)), not the column radius $r_in"))
    dz = domain[:dx]
    ring(c) = (j = cld(c, nz); π * (radii[j + 1]^2 - radii[j]^2))
    domain[:volumes] = [ring(c) * dz[c] for c in 1:(nz * nr)]
    N, n, fc = domain[:neighbors], domain[:normals], domain[:face_centroids]
    domain[:areas, Jutul.Faces()] = [abs(n[1, f]) > 0.5 ? ring(N[1, f]) : 2π * fc[2, f] * dz[N[1, f]] for f in axes(N, 2)]
    bn, bc, bnb = domain[:boundary_normals], domain[:boundary_centroids], domain[:boundary_neighbors]
    domain[:boundary_areas, Jutul.BoundaryFaces()] = [abs(bn[1, f]) > 0.5 ? ring(bnb[f]) :
        (abs(bn[2, f]) > 0.5 ? 2π * bc[2, f] * dz[bnb[f]] : 0.0) for f in eachindex(bnb)]
    domain[:radial_cells, Unit()] = nr
    return domain
end

# Whether a data domain is the (z, r) plane of a column, set up by mocca_domain
# with axisymmetric = true
_is_axisymmetric(data_domain) = haskey(data_domain, :radial_cells, Unit())

# Whether each face lies across the axis, between two cells along it
_axial_faces(data_domain) = abs.(data_domain[:normals][1, :]) .> 0.5

# The cells the wall touches: the outer ring of an axisymmetric column, and
# every cell otherwise
function _wall_cell_mask(data_domain)
    nc = Jutul.number_of_cells(Jutul.physical_representation(data_domain))
    _is_axisymmetric(data_domain) || return trues(nc)
    nr = Int(first(data_domain[:radial_cells, Unit()]))
    nz = nc ÷ nr
    return [cld(c, nz) == nr for c in 1:nc]
end

"""
    wall_cells(model) → Vector{Int}

The cells of a distributed unit that its wall touches: the outer ring of a
column meshed across its radius, and every cell of one that is not.
"""
wall_cells(model) = findall(_wall_cell_mask(model.data_domain))

# ------------------------------------------------------------------------------

# Inner or outer surface [m²] of the wall around each cell; zero for the cells
# inside an axisymmetric column, which the wall does not touch
struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    dx = data_domain[:dx]
    if T == :in
        r = first(data_domain[:r_in, Unit()])
    else
        @assert T == :out
        r = first(data_domain[:r_out, Unit()])
    end
    return 2π .* r .* dx .* _wall_cell_mask(data_domain)  # Lateral surface area of cylindrical shell
end

struct CellDx <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::CellDx, symb)
    return copy(data_domain[:dx])
end

# Cell volumes (solid + pore)
struct BulkVolume <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::BulkVolume, symb)
    if haskey(data_domain, :volumes)
        return copy(data_domain[:volumes])
    else
        error(":volumes must be present in DataDomain to initialize $symb")
    end
end

# Pore volume = volume × porosity
struct FluidVolume <: Jutul.ScalarVariable end
Jutul.minimum_value(::FluidVolume) = eps()

function Jutul.default_parameter_values(data_domain, model, param::FluidVolume, symb)
    if haskey(data_domain, :fluid_volume, Jutul.Cells())
        return copy(data_domain[:fluid_volume])
    elseif haskey(data_domain, :volumes) && haskey(data_domain, :porosity)
        return copy(data_domain[:volumes] .* data_domain[:porosity])
    else
        error(":fluid_volume or (:volumes + :porosity) must be present in DataDomain to initialize $symb")
    end
end

# Cell permeability [m²], used for the half-cell connection at a port
struct Permeability <: Jutul.ScalarVariable end
Jutul.variable_scale(::Permeability) = 1e-10
Jutul.minimum_value(::Permeability) = 0.0

function Jutul.default_parameter_values(data_domain, model, param::Permeability, symb)
    return copy(data_domain[:permeability])
end

# Flow cross-section of the bed [m²]
struct BedCrossSectionArea <: Jutul.ScalarVariable end
Jutul.associated_entity(::BedCrossSectionArea) = Unit()
Jutul.minimum_value(::BedCrossSectionArea) = 0.0

function Jutul.default_parameter_values(data_domain, model, ::BedCrossSectionArea, symb)
    r_in = first(data_domain[:r_in, Unit()])
    return [π * r_in^2]
end

"""
    compute_column_face_area(model, state)

Flow cross-section of a distributed unit [m²]. Read from the `BedCrossSectionArea`
parameter when `state` carries parameters, otherwise from the column radius in
the data domain.
"""
function compute_column_face_area(model::DistributedUnitModel, state)
    if _state_has(state, :BedCrossSectionArea)
        return state.BedCrossSectionArea[1]
    else
        r_in = first(model.data_domain[:r_in, Unit()])
        return π * r_in^2
    end
end

"""
    calc_bc_trans(model, state, cell)

Transmissibility [m³] of the half-cell connection between `cell` and the end of the
unit it touches.
"""
function calc_bc_trans(model::DistributedUnitModel, state, cell)
    k = _state_or_domain(state, model, :Permeability, :permeability)[cell]
    dx = _state_or_domain(state, model, :CellDx, :dx)[cell] / 2.0
    A = compute_column_face_area(model, state)
    return k * A / dx
end

# TPFA transmissibilities (face parameter)
struct Transmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::Transmissibilities) = 1e-10
Jutul.minimum_value(::Transmissibilities) = 0.0
Jutul.associated_entity(::Transmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::Transmissibilities, symb)
    if haskey(data_domain, :transmissibilities, Jutul.Faces())
        return copy(data_domain[:transmissibilities])
    elseif haskey(data_domain, :permeability, Jutul.Cells())
        return Jutul.compute_face_trans(data_domain, data_domain[:permeability])
    else
        error(":permeability or :transmissibilities must be present in DataDomain to initialize $symb")
    end
end

# Inject Transmissibilities when PotentialFlow discretization is present
function Jutul.select_parameters!(S, disc::Jutul.PotentialFlow, model::Jutul.SimulationModel{D, S2}) where {D, S2 <: DistributedUnit}
    S[:Transmissibilities] = Transmissibilities()
end
