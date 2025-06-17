export mocca_domain
"""
    mocca_domain(mesh::Jutul.CartesianMesh, system::AdsorptionSystem)

Set up a data domain for carbon capture simulation.
"""
function mocca_domain(mesh, system::AdsorptionSystem; kwarg...)
    domain = JutulDarcy.reservoir_domain(mesh, porosity = system.p.Φ, permeability = system.permeability)
    domain[:diffusion_coefficient] = system.dispersion
    domain[:thermal_conductivity] = system.p.K_z
    nc = Jutul.number_of_cells(mesh)

    dx = map(i -> first(Jutul.cell_dims(mesh, i)), 1:nc)
    domain[:dx] = dx

    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end

# TODO: This causes problem for adjoint simulation. Is it needed?
function Jutul.select_linear_solver(model::AdsorptionModel; kwarg...)
    #return Jutul.LUSolver(; kwarg...)
    return nothing
end
