"""
    Unit <: Jutul.JutulEntity

Entity holding one set of scalar values per unit operation (bed radius, wall
properties, fluid viscosity, ...). Every Mocca domain registers it with count 1,
and the matching Jutul parameters are read at runtime as `state.ParamName[1]`.

`Column` is kept as an alias for backwards compatibility.
"""
struct Unit <: Jutul.JutulEntity end
const Column = Unit

"""
    MoccaSystem <: Jutul.JutulSystem

Common supertype for every Mocca unit operation. Each unit in a flowsheet is a
`Jutul.SimulationModel` whose system is a subtype of one of:

- [`DistributedUnit`](@ref): has a mesh and solves PDEs (beds, membrane sides, ...)
- [`LumpedUnit`](@ref): 0D with holdup (tanks, nodes, flash drums)
- [`FlowElement`](@ref): no holdup, owns a flow and solves a device law (valves,
  flow controllers, pumps, compressors)

Implementations must provide a `component_names` field or override
[`component_names`](@ref).
"""
abstract type MoccaSystem <: Jutul.JutulSystem end

"Unit with a mesh that solves conservation laws on cells and faces."
abstract type DistributedUnit <: MoccaSystem end

"Distributed unit where gas and a sorbent share each cell."
abstract type SorbentBed <: DistributedUnit end

"Distributed unit forming one side of a two-channel device (membrane, absorber, heat exchanger)."
abstract type FlowChannel <: DistributedUnit end

"0D unit with holdup."
abstract type LumpedUnit <: MoccaSystem end

"Unit without holdup that owns a flow and solves a device law."
abstract type FlowElement <: MoccaSystem end

const MoccaModel = Jutul.SimulationModel{<:Any, <:MoccaSystem, <:Any, <:Any}
const DistributedUnitModel = Jutul.SimulationModel{<:Any, <:DistributedUnit, <:Any, <:Any}
const SorbentBedModel = Jutul.SimulationModel{<:Any, <:SorbentBed, <:Any, <:Any}
const FlowElementModel = Jutul.SimulationModel{<:Any, <:FlowElement, <:Any, <:Any}

"""
    component_names(sys::MoccaSystem)

Names of the chemical components tracked by `sys`, in index order.
"""
component_names(sys::MoccaSystem) = sys.component_names

"""
    number_of_components(sys::MoccaSystem)

Number of chemical components tracked by `sys`.
"""
number_of_components(sys::MoccaSystem) = length(component_names(sys))

# Set up a discretized domain with TPFA potential flow for distributed units.
function Jutul.discretize_domain(d::Jutul.DataDomain, system::DistributedUnit, ::Val{:default}; general_ad = true, kwarg...)
    g = Jutul.physical_representation(d)
    N = d[:neighbors]
    nc = Jutul.number_of_cells(g)
    flow_disc = Jutul.PotentialFlow(N, nc)
    disc = (mass_flow = flow_disc,)
    domain = Jutul.DiscretizedDomain(g, disc; kwarg...)
    Jutul.transfer_entities!(domain, d)
    return domain
end

# Mocca problems are small, so a direct sparse solve is used throughout.
# TODO: This causes problem for adjoint simulation. Is it needed?
function Jutul.select_linear_solver(model::MoccaModel; kwarg...)
    return nothing
end

# Read a value that is a Jutul parameter during simulation, but may be missing
# from a state passed to post-processing (e.g. an objective function). Falls
# back to the model's data domain.
function _state_or_domain(state, model, symb::Symbol, key::Symbol, entity = Jutul.Cells())
    if _state_has(state, symb)
        return _state_get(state, symb)
    else
        return model.data_domain[key, entity]
    end
end

_state_has(state::AbstractDict, k) = haskey(state, k)
_state_has(state, k) = hasproperty(state, k)
_state_get(state::AbstractDict, k) = state[k]
_state_get(state, k) = getproperty(state, k)
