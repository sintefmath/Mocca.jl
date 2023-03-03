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

compute_permeability(sys::AdsorptionFlowSystem) = 4 / 150 * ((sys.p.Φ / (1 - sys.p.Φ))^2) * (sys.p.d_p / 2)^2 * sys.p.Φ
axial_dispersion(sys::AdsorptionFlowSystem) = 0.7 * sys.p.D_m + 0.5 * sys.p.V0_inter * sys.p.Φ * sys.p.d_p

"Area of column wall [m^2]"
area_wall(sys::AdsorptionFlowSystem) = π * (sys.p.r_out^2 - sys.p.r_in^2)