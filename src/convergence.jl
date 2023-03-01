using LinearAlgebra
using Tullio
function Jutul.convergence_criterion(model::Jutul.SimulationModel{D,S}, storage, eq::Jutul.ConservationLaw, eq_s, r; dt=1) where {D,S<:AdsorptionFlowSystem}
    M = Jutul.global_map(model.domain)

    v = x -> Jutul.as_value(Jutul.active_view(x, M, for_variables=false))
    #Φ = v(storage.state.FluidVolume)
    Φ = model.system.Φ

    scale = 1.0
    # scale = 1/1000.0
    # scale = 1/1e13
    e = scale .* [maximum(abs.(r[j, :]) * dt / (Φ)) for j in 1:2]

    # e = sum(abs, r, dims = 2)
    y = storage.state.y
    q = storage.state.adsorptionRates
    # @show r
    # @info "conv " size(r)
    #r_scaled = [r[:, cell] ./ [y[:, cell]..., q[:, cell]...] for cell in eachindex(q[1, :])]

    e = norm(r)

    N = length(Φ)
    pv_t = sum(Φ)
    avg_density = 1.0#sum(ρ, dims = 2)./N
    r_sum = sum(r, dims=2)
    # mb = @. (dt/pv_t)*abs(r_sum)/avg_density

    names = ["CO2", "N2"]
    R = (
        CNV=(errors=e, names=names),
        # MB = (errors = mb, names = names)
    )
    return R
end