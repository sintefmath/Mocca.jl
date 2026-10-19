# Packed-bed correlations for the parameters of a fixed bed.

function compute_permeability(porosity, d_p)
    return 4 / 150 * ((porosity / (1 - porosity))^2) * (d_p / 2)^2 * porosity
end

function compute_dispersion(D_m, V0_inter, d_p)
    return 0.7 * D_m + 0.5 * V0_inter * d_p
end
