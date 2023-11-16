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
    dx = zeros(nc)
    for i in 1:nc
        dx[i] = first(Jutul.cell_dims(mesh, i))
    end
    domain[:dx] = dx

    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end
