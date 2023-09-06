
abstract type AdsorptionVelocityModel end


@with_kw struct PlugFlow{Φ, d_p} <: AdsorptionVelocityModel
    
    permeability = 4 / 150 * ((Φ / (1 - Φ))^2) * (d_p / 2)^2 * Φ

end

abstract type AdsorptionDispersionModel end

@with_kw struct StandardDispersion{D_m, V0_inter, d_p} <: AdsorptionDispersionModel
    
    dispersion = 0.7 * p.D_m + 0.5 * p.V0_inter * p.d_p

end