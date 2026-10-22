# Building flowsheets from materials: beds from an adsorbent and a column, and
# devices with roles.

"""
How a bed was built, so that it can be rebuilt and written to a file.
"""
mutable struct BedSpec
    adsorbent::Adsorbent
    column::Column
    ncells::Int
    thermal::AbstractThermalModel
    dispersion_velocity::Union{Symbol, Float64}
    dispersion::Union{Nothing, Float64}
    permeability::Union{Nothing, Float64}
    velocity_used::Float64
    radial_cells::Int
    radii::Union{Nothing, Vector{Float64}}
    radial_dispersion::Union{Nothing, Float64}
end

"""
    Flowsheet(gas::Gas, feed::Feed)

An empty flowsheet for `gas`, fed with `feed`. Add beds with
[`add_bed!`](@ref) and devices with [`add_device!`](@ref), and run it as a
[`Process`](@ref).

Devices have roles (feed, product, extract, ...) that give their default
boundary conditions and valve sizes, tell the step templates which device to
open, and tell [`performance`](@ref) which streams are feed and product. The
units and connections are a [`UnitNetwork`](@ref), so the lower-level
functions, such as [`add_unit!`](@ref), [`connect!`](@ref) and
[`setup_flowsheet_model`](@ref), work on a flowsheet too.
"""
struct Flowsheet <: AbstractFlowsheet
    network::UnitNetwork
    gas::Gas
    feed::Feed
    roles::Dict{Symbol, Symbol}            # device => role
    specs::Dict{Symbol, BedSpec}           # bed => how it was built
    default_boundaries::Set{PortRef}       # boundaries set from a role's defaults
end
Flowsheet(gas::Gas, feed::Feed) = Flowsheet(UnitNetwork(), gas, feed, Dict{Symbol, Symbol}(), Dict{Symbol, BedSpec}(), Set{PortRef}())

Networks.network(fs::Flowsheet) = fs.network

function Networks.set_boundary!(fs::Flowsheet, device::Symbol, port::Symbol, condition::PortCondition)
    set_boundary!(fs.network, device, port, condition)
    delete!(fs.default_boundaries, PortRef(device, port))
    return fs
end

# ------------------------------------------------------------------------------
# Beds
# ------------------------------------------------------------------------------

"""
    fixed_bed_model(gas, feed, adsorbent, column; ncells = 30, thermal = WithWall(),
        dispersion_velocity = 1.0, dispersion = nothing, permeability = nothing,
        radial_cells = 1, radii = nothing, radial_dispersion = nothing)

The model of a fixed bed of `adsorbent` packed in `column`, for `gas`, on
`ncells` cells along the column. The permeability (Kozeny–Carman, from the
bed porosity and particle size) and the axial dispersion (`0.7 D_m + 0.5 v
d_p`, at the interstitial `dispersion_velocity` [m/s]) are computed unless
given. The gas viscosity defaults to that of the `feed`; see [`Gas`](@ref).

With `radial_cells` above 1, each of the cells along the column is split into
rings, uniform in radius unless `radii` gives their edges from 0 to the
column radius, and the gas disperses and conducts across the radius as well
(see [`column_mesh`](@ref) and [`mocca_domain`](@ref)). The radial dispersion
is `0.7 D_m + 0.1 v d_p` (see [`compute_radial_dispersion`](@ref)) unless
`radial_dispersion` gives it, and the column must give its
`radial_conductivity` and, for a bed with a wall, its `wall_heat_transfer`,
which takes the place of `inner_heat_transfer` (see [`Column`](@ref)).
"""
function fixed_bed_model(gas::Gas, feed::Feed, adsorbent::Adsorbent, column::Column; ncells = 30, thermal = WithWall(),
        dispersion_velocity = 1.0, dispersion = nothing, permeability = nothing,
        radial_cells = 1, radii = nothing, radial_dispersion = nothing)
    loc = isempty(adsorbent.name) ? "bed" : "bed of $(adsorbent.name)"
    ncells >= 2 || input_error(loc, "ncells", "needs at least 2 cells, got $ncells")
    radial_cells >= 1 || input_error(loc, "radial_cells", "needs at least 1 cell, got $radial_cells")
    thermal isa WithWall && column.outer_radius <= column.radius &&
        input_error(loc, "column", "a bed with a wall needs the column's outer_radius; use thermal = Adiabatic() for a bed without one")
    radial = radial_cells > 1
    if radial
        has_energy_balance_ = !(thermal isa Isothermal)
        has_energy_balance_ && isnothing(column.radial_conductivity) &&
            input_error(loc, "column", "a bed with radial_cells needs the column's radial_conductivity")
        thermal isa WithWall && isnothing(column.wall_heat_transfer) &&
            input_error(loc, "column", "a bed with radial_cells and a wall needs the column's wall_heat_transfer")
        isnothing(radii) || length(radii) == radial_cells + 1 ||
            input_error(loc, "radii", "needs radial_cells + 1 = $(radial_cells + 1) ring edges, got $(length(radii))")
    else
        isnothing(radii) || input_error(loc, "radii", "are only used with radial_cells above 1")
    end
    system = FixedBed(;
        isotherm = isotherm(adsorbent, gas),
        mass_transfer = mass_transfer(adsorbent, gas),
        molecular_masses = gas.molar_masses,
        component_names = gas.names,
        heat_capacity_gas = gas.heat_capacities,
        heat_capacity_adsorbed = adsorbed_heat_capacities(adsorbent, gas),
        thermal = thermal)
    perm = isnothing(permeability) ? compute_permeability(column.porosity, adsorbent.particle_diameter) : permeability
    disp = isnothing(dispersion) ? compute_dispersion(axial_diffusivity(gas), dispersion_velocity, adsorbent.particle_diameter) : dispersion
    mesh = radial ? column_mesh(ncells, column.length, column.radius; radial_cells = radial_cells, radii = radii) :
        column_mesh(ncells, column.length, column.radius)
    radial_kw = radial ? (axisymmetric = true,
        radial_diffusion_coefficient = isnothing(radial_dispersion) ?
            compute_radial_dispersion(axial_diffusivity(gas), dispersion_velocity, adsorbent.particle_diameter) : radial_dispersion,
        radial_thermal_conductivity = something(column.radial_conductivity, column.axial_conductivity)) : (;)
    h_in = radial && thermal isa WithWall ? column.wall_heat_transfer : column.inner_heat_transfer
    domain = mocca_domain(mesh; radial_kw...,
        porosity = column.porosity, permeability = perm,
        r_in = column.radius, r_out = column.outer_radius,
        diffusion_coefficient = disp, thermal_conductivity = column.axial_conductivity,
        adsorbent_density = adsorbent.density,
        adsorbent_heat_capacity = adsorbent.heat_capacity,
        wall_density = column.wall.density,
        wall_heat_capacity = column.wall.heat_capacity,
        wall_conductivity = column.wall.conductivity,
        fluid_viscosity = gas_viscosity(gas, feed.composition),
        inner_htc = h_in,
        outer_htc = column.outer_heat_transfer,
        ambient_temperature = column.ambient_temperature)
    return setup_process_model(system, domain)
end

"""
    add_bed!(fs, name, adsorbent, column; ncells = 30, thermal = WithWall(),
        dispersion_velocity = :adsorption, dispersion = nothing, permeability = nothing,
        radial_cells = 1, radii = nothing, radial_dispersion = nothing)

Add a fixed bed of `adsorbent` in `column` to the flowsheet under `name`. Its
ports are `:bottom` and `:top`.

The axial dispersion depends on the gas velocity. By default
(`dispersion_velocity = :adsorption`) it is set from the interstitial velocity
of the bed's adsorption step when a [`Process`](@ref) is built, or the largest
adsorption velocity in the cycle for a bed fed from another; it then stays
fixed as settings change. Give a velocity [m/s], or the `dispersion` [m²/s]
itself, to fix it. It is the bed's `AxialDispersion` parameter, so it can also
be set through `parameters` of a `Process` or a design variable.

With `radial_cells` above 1 the bed is meshed across its radius too, in rings
from the axis to the wall, and resolves the gradients across it: the column
must then give its `radial_conductivity` and `wall_heat_transfer` (see
[`Column`](@ref)). Its ports are the end rings, so the devices, steps and
results are those of a bed of one radial cell. The radial dispersion follows
the same velocity as the axial one, unless `radial_dispersion` [m²/s] fixes
it, and is the bed's `RadialDispersion` parameter. See
[`fixed_bed_model`](@ref) for the other keywords.
"""
function add_bed!(fs::Flowsheet, name::Symbol, adsorbent::Adsorbent, column::Column; ncells = 30, thermal = WithWall(),
        dispersion_velocity = :adsorption, dispersion = nothing, permeability = nothing,
        radial_cells = 1, radii = nothing, radial_dispersion = nothing)
    dispersion_velocity isa Symbol && dispersion_velocity != :adsorption &&
        input_error("bed $name", "dispersion_velocity", "use :adsorption or a velocity in m/s, got $dispersion_velocity")
    v = dispersion_velocity isa Symbol ? 1.0 : Float64(dispersion_velocity)
    model = _at(() -> fixed_bed_model(fs.gas, fs.feed, adsorbent, column; ncells = ncells, thermal = thermal,
        dispersion_velocity = v, dispersion = dispersion, permeability = permeability,
        radial_cells = radial_cells, radii = radii, radial_dispersion = radial_dispersion), "bed $name")
    add_unit!(fs, name, model)
    fs.specs[name] = BedSpec(adsorbent, column, ncells, thermal,
        dispersion_velocity isa Symbol ? dispersion_velocity : Float64(dispersion_velocity),
        isnothing(dispersion) ? nothing : Float64(dispersion), isnothing(permeability) ? nothing : Float64(permeability), v,
        radial_cells, isnothing(radii) ? nothing : collect(Float64, radii),
        isnothing(radial_dispersion) ? nothing : Float64(radial_dispersion))
    return fs
end

# The bed as add_bed! built it from `spec`, on `ncells` cells along the column
# and `radial_cells` across it
function _bed_model(fs::Flowsheet, spec::BedSpec; ncells = spec.ncells, radial_cells = spec.radial_cells)
    radii = radial_cells == spec.radial_cells ? spec.radii : nothing
    return fixed_bed_model(fs.gas, fs.feed, spec.adsorbent, spec.column; ncells = ncells, thermal = spec.thermal,
        dispersion_velocity = spec.velocity_used, dispersion = spec.dispersion, permeability = spec.permeability,
        radial_cells = radial_cells, radii = radii, radial_dispersion = spec.radial_dispersion)
end

# The cells of a bed, as in "30" or "30 × 4" for one meshed across its radius
_cells_of(spec::BedSpec) = spec.radial_cells > 1 ? "$(spec.ncells) × $(spec.radial_cells)" : "$(spec.ncells)"

"""
    bed_dispersion(fs, name) → Float64

Axial dispersion [m²/s] of bed `name`, added with [`add_bed!`](@ref): the
`dispersion` given, or `0.7 D_m + 0.5 v d_p` at the velocity it was set for.
"""
function bed_dispersion(fs::Flowsheet, name::Symbol)
    spec = fs.specs[name]
    isnothing(spec.dispersion) || return spec.dispersion
    return compute_dispersion(axial_diffusivity(fs.gas), spec.velocity_used, spec.adsorbent.particle_diameter)
end

"""
    bed_radial_dispersion(fs, name) → Float64

Radial dispersion [m²/s] of bed `name`, added with [`add_bed!`](@ref) with
`radial_cells` above 1: the `radial_dispersion` given, or
`0.7 D_m + 0.1 v d_p` at the velocity the axial dispersion was set for.
"""
function bed_radial_dispersion(fs::Flowsheet, name::Symbol)
    spec = fs.specs[name]
    isnothing(spec.radial_dispersion) || return spec.radial_dispersion
    return compute_radial_dispersion(axial_diffusivity(fs.gas), spec.velocity_used, spec.adsorbent.particle_diameter)
end

"Units of `fs` that hold sorbent (see [`holds_sorbent`](@ref)), in flowsheet order."
beds(fs::AbstractFlowsheet) = [n for n in unit_names(fs) if holds_sorbent(fs[n])]

# ------------------------------------------------------------------------------
# Devices with roles
# ------------------------------------------------------------------------------

"""
Roles of flow devices:

| Role           | Connected by | Boundary | Meaning |
|----------------|--------------|----------|---------|
| `:feed`        | its outlet   | inlet    | feed gas into a bed; counts as feed |
| `:light_inlet` | its outlet   | inlet    | light product into a bed, for light product pressurisation or purge |
| `:product`     | its inlet    | outlet   | light product (raffinate) out of a bed |
| `:extract`     | its inlet    | outlet   | heavy product out of a bed, through a vacuum pump; counts as product |
| `:waste`       | its inlet    | outlet   | other gas out of a bed, such as a water purge |
| `:link`        | both ports   | none     | between two beds: series feed, equalisation |
"""
const DEVICE_ROLES = (:feed, :light_inlet, :product, :extract, :waste, :link)
const INLET_ROLES = (:feed, :light_inlet)
const OUTLET_ROLES = (:product, :extract, :waste)

_port_pair(x::Tuple{Symbol, Symbol}) = x
_port_pair(x::Pair{Symbol, Symbol}) = (first(x), last(x))
_port_pair(x) = input_error("device", "", "a connection is (unit, port), such as (:Bed, :bottom); got $x")

"""
    add_device!(fs, name, role; at, boundary = nothing)
    add_device!(fs, name, :link; from, to)

Add a flow device with a `role` (see [`DEVICE_ROLES`](@ref)), connected to a
bed port `at = (bed, port)`, such as `(:Bed, :bottom)`. A link joins `from`,
a bed port, to `to`, another.

The boundary of an inlet or outlet device defaults to the feed temperature, a
role composition (the feed for feed, extract and waste; the light product for
product and light inlet) and a pressure that steps override: the feed
pressure, or for extract and waste the lowest pressure of the bed in its
cycle. Pass `boundary` (a [`PortCondition`](@ref), or a NamedTuple with any of
`pressure`, `temperature` and `composition`) to set it.
"""
function add_device!(fs::Flowsheet, name::Symbol, role::Symbol; at = nothing, from = nothing, to = nothing, boundary = nothing)
    loc = "device $name"
    role in DEVICE_ROLES || input_error(loc, "role", "unknown role $role; roles are $(join(DEVICE_ROLES, ", "))")
    add_unit!(fs, name, setup_device_model(fs.gas.names))
    fs.roles[name] = role
    try
        _at(loc) do
            if role == :link
                (isnothing(from) || isnothing(to)) && input_error(loc, "", "a link needs from = (bed, port) and to = (bed, port)")
                isnothing(boundary) || input_error(loc, "boundary", "a link has no boundary")
                connect!(fs, _port_pair(from) => (name, :inlet))
                connect!(fs, _port_pair(to) => (name, :outlet))
            else
                isnothing(at) && input_error(loc, "at", "give the bed port it connects to, such as at = (:Bed, :bottom)")
                inlet = role in INLET_ROLES
                connect!(fs, _port_pair(at) => (name, inlet ? :outlet : :inlet))
                port = inlet ? :inlet : :outlet
                set_boundary!(fs.network, name, port, _role_boundary(fs, role, boundary, loc))
                isnothing(boundary) && push!(fs.default_boundaries, PortRef(name, port))
            end
        end
    catch
        # Leave the flowsheet as it was, so that the call can be corrected
        # and made again under the same name
        _remove_device!(fs, name)
        rethrow()
    end
    return fs
end

# Undo add_device!: the device, its role, connections and boundaries
function _remove_device!(fs::Flowsheet, name::Symbol)
    net = fs.network
    delete!(net.units, name)
    filter!(!=(name), net.order)
    filter!(c -> all(r -> r.name != name, connection_ports(c)), net.connections)
    for p in DEVICE_PORTS
        delete!(net.boundaries, PortRef(name, p))
        delete!(fs.default_boundaries, PortRef(name, p))
    end
    delete!(fs.roles, name)
    return fs
end

role_composition(fs::Flowsheet, role::Symbol) = role in (:product, :light_inlet) ? fs.feed.light_product : fs.feed.composition

function _role_boundary(fs, role, boundary, loc; pressure = fs.feed.pressure)
    default = PortCondition(pressure = pressure, temperature = fs.feed.temperature, composition = role_composition(fs, role))
    isnothing(boundary) && return default
    boundary isa PortCondition && return boundary
    p = get(boundary, :pressure, default.pressure)
    T = get(boundary, :temperature, default.temperature)
    y = haskey(boundary, :composition) ? composition(fs.gas, boundary.composition; location = loc) : default.composition
    return PortCondition(pressure = p, temperature = T, composition = y)
end

"Role of device `d` (see [`add_device!`](@ref)), or `nothing`."
device_role(fs::AbstractFlowsheet, d::Symbol) = nothing
device_role(fs::Flowsheet, d::Symbol) = get(fs.roles, d, nothing)
devices_with_role(fs::AbstractFlowsheet, role::Symbol) = [d for d in unit_names(fs) if device_role(fs, d) == role]

"Devices connected at port `port` of unit `u`, in flowsheet order."
function devices_at(fs::AbstractFlowsheet, u::Symbol, port::Symbol)
    r = PortRef(u, port)
    return [d for d in unit_names(fs) if any(c -> c.unit == r && c.device.name == d, stream_connections(fs))]
end

"The [`StreamConnection`](@ref)s of device `d`."
attachments(fs::AbstractFlowsheet, d::Symbol) = [c for c in stream_connections(fs) if c.device.name == d]

"""
    default_conductance(fs, device)

Valve conductance [m³/(s·Pa)] for `device`: that of the half cell at the bed
end it is attached to (see [`port_conductance`](@ref)), so that the valve adds
no more resistance than the bed's end and the bed follows the pressure set at
the other side. For a link, the two half cells in series.
"""
function default_conductance(fs::Flowsheet, d::Symbol)
    C = [port_conductance(fs[c.unit.name], c.unit.port) for c in attachments(fs, d)]
    isempty(C) && input_error("device $d", "", "it is not connected to any unit")
    return length(C) == 1 ? C[1] : 1 / (1 / C[1] + 1 / C[2])
end

"Volumetric feed rate [m³/s] into `bed` for an interstitial velocity `v` [m/s]."
function interstitial_flow_rate(bed, v)
    A = compute_column_face_area(bed, NamedTuple(Jutul.setup_parameters(bed)))
    ε = first(bed.data_domain[:porosity])
    return ε * v * A
end

"Volumetric feed rate [m³/s] into `bed` for a superficial velocity `u` [m/s]."
superficial_flow_rate(bed, u) = u * compute_column_face_area(bed, NamedTuple(Jutul.setup_parameters(bed)))

# ------------------------------------------------------------------------------
# Inspection
# ------------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", fs::AbstractFlowsheet)
    nb = length(holdup_units(fs))
    nd = length(unit_names(fs)) - nb
    println(io, nameof(typeof(fs)), " with $nb unit$(nb == 1 ? "" : "s") and $nd device$(nd == 1 ? "" : "s")",
        fs isa Flowsheet ? " for $(join(fs.gas.names, ", "))" : "")
    print_flowsheet_diagram(io, fs)
end

_unit_spec(fs::AbstractFlowsheet, u) = nothing
_unit_spec(fs::Flowsheet, u) = get(fs.specs, u, nothing)

"""
    print_flowsheet_diagram([io], fs)

Print each unit of `fs` with the devices at each of its ports, their roles and
boundaries.
"""
print_flowsheet_diagram(fs::AbstractFlowsheet) = print_flowsheet_diagram(stdout, fs)
function print_flowsheet_diagram(io::IO, fs::AbstractFlowsheet)
    for u in holdup_units(fs)
        spec = _unit_spec(fs, u)
        desc = spec isa BedSpec ? "fixed bed of $(isempty(spec.adsorbent.name) ? "adsorbent" : spec.adsorbent.name), $(_cells_of(spec)) cells" : string(nameof(typeof(fs[u].system)))
        println(io, "  ", u, ": ", desc)
        for p in keys(ports(fs[u]))
            ds = devices_at(fs, u, p)
            isempty(ds) && continue
            println(io, "    ", rpad(string(p), 7), join([_device_label(fs, d, u) for d in ds], "; "))
        end
    end
end

function _device_label(fs, d, from_unit)
    role = device_role(fs, d)
    s = string(d, isnothing(role) ? "" : " ($role)")
    for c in attachments(fs, d)
        c.unit.name == from_unit || (s *= " to $(c.unit.name) $(c.unit.port)")
    end
    for p in DEVICE_PORTS
        haskey(network(fs).boundaries, PortRef(d, p)) || continue
        c = boundary_condition(fs, d, p)
        s *= ", $(p == :inlet ? "from" : "to") boundary at $(round(evaluate_profile(c.pressure, 0.0) / 1e5, sigdigits = 4)) bar"
    end
    return s
end
