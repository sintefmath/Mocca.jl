export mocca_domain
"""
    mocca_domain(mesh::Jutul.CartesianMesh, system::AdsorptionSystem)

Set up a data domain for carbon capture simulation. Currently only implemented
for a CartesianMesh instance.
"""
function mocca_domain(mesh::Jutul.CartesianMesh, system::AdsorptionSystem; kwarg...)
    domain = JutulDarcy.reservoir_domain(mesh, porosity = system.p.Φ, permeability = system.permeability)
    domain[:diffusion_coefficient] = system.dispersion
    domain[:thermal_conductivity] = system.p.K_z 
    domain[:dx] = mesh.dims[2]
    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end
