using Parameters
@with_kw struct AdsorptionFlowSystem <: JutulDarcy.MultiComponentSystem
    """
    Number of components of the system. For now we only support two (CO2 and N2).

    Some rewriting is expected when changing this. 
    """
    number_of_components::Int64 = 2
    
    """
    Use this to turn on and off the forcing term in the equation.
    If you don't know what this does, leave it to 1.0.
    """
    forcing_term_coefficient::Float64 = 1.0

    "Contains all the relevant physical constants for the system"
    p::AdsorptionParameters = AdsorptionParameters()
end
const AdsorptionFlowModel = Jutul.SimulationModel{<:Any,<:AdsorptionFlowSystem,<:Any,<:Any}

JutulDarcy.number_of_components(sys::AdsorptionFlowSystem) = sys.number_of_components
JutulDarcy.has_other_phase(::AdsorptionFlowSystem) = false



compute_permeability(sys::AdsorptionFlowSystem) = compute_permeability(sys.p)
axial_dispersion(sys::AdsorptionFlowSystem) = axial_dispersion(sys.p)

function compute_dx(model::AdsorptionFlowModel, self_cell)\
    # TODO: We need to get dx in a nicer way
    g = JutulDarcy.physical_representation(model.data_domain)

    return first(g.deltas)
end
"Area of column wall [m^2]"
area_wall(sys::AdsorptionFlowSystem) = π * (sys.p.r_out^2 - sys.p.r_in^2)

area_wall_in(sys::AdsorptionFlowSystem, Δx) = π * sys.p.r_in * 2 * Δx
area_wall_in(model::AdsorptionFlowModel, self_cell) = area_wall_in(model.system, compute_dx(model, self_cell))

area_wall_out(sys::AdsorptionFlowSystem, Δx) = π * sys.p.r_out * 2 * Δx
area_wall_out(model::AdsorptionFlowModel, self_cell) = area_wall_out(model.system, compute_dx(model, self_cell))