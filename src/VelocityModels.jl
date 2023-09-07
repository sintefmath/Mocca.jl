#TODO: check this is the best way for interchangeable velocity models, dispersion models

function calc_perm_plug_flow(Φ::Float64, d_p::Float64) 
    return 4 / 150 * ((Φ / (1 - Φ))^2) * (d_p / 2)^2 * Φ
end
function calc_dispersion(D_m::Float64, V0_inter::Float64, d_p::Float64)
    return 0.7 * D_m + 0.5 * V0_inter * d_p
end