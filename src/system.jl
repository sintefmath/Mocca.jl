# using Parameters
# @with_kw struct AdsorptionSystem <: JutulDarcy.MultiComponentSystem
#     """
#     Number of components of the system. For now we only support two (CO2 and N2).

#     Some rewriting is expected when changing this. 
#     """
#     number_of_components::Int64 = 2
    
#     """
#     Use this to turn on and off the forcing term in the equation.
#     If you don't know what this does, leave it to 1.0.
#     """
#     forcing_term_coefficient::Float64 = 1.0

#     "Contains all the relevant physical constants for the system"
#     p::HaghpanahParameters = HaghpanahParameters() #TODO: is this necessary?
# end
