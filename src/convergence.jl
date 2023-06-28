using LinearAlgebra
using Tullio
function Jutul.convergence_criterion(model::Jutul.SimulationModel{D,S}, storage, eq::Jutul.ConservationLaw, eq_s, r; dt=1) where {D,S<:AdsorptionFlowSystem}

    # # TODO: Clean this up and fix it. 
    # M = Jutul.global_map(model.domain)

    # v = x -> Jutul.as_value(Jutul.active_view(x, M, for_variables=false))
    # #Φ = v(storage.state.FluidVolume)
    # Φ = model.system.p.Φ

    # scale = 1.0
    # # scale = 1/1000.0
    # # scale = 1/1e13
    # #e = scale .* [maximum(abs.(r[j, :]) * dt / (Φ)) for j in 1:2]

    # # e = sum(abs, r, dims = 2)
    # y = storage.state.y
    # q = storage.state.AdsorbedConcentration
    # # @show r
    # # @info "conv " size(r)
    # #r_scaled = [r[:, cell] ./ [y[:, cell]..., q[:, cell]...] for cell in eachindex(q[1, :])]

    # e = norm(r, Inf)#/1e5

    # N = length(Φ)
    # pv_t = sum(Φ)
    # avg_density = 1.0#sum(ρ, dims = 2)./N
    # r_sum = sum(r, dims=2)
    # # mb = @. (dt/pv_t)*abs(r_sum)/avg_density

    # names = ["CO2", "N2"]
    # n = size(r, 1)
    # e = zeros(2)
    # for i in eachindex(names)
    #     e[i] = norm(r[i, :], Inf)
    # end
    # R = (
    #     CNV=(errors=e, names=names),
    #     # MB = (errors = mb, names = names)
    # )
    n = Jutul.number_of_equations_per_entity(model, eq)
    @tullio max e[i] := abs(r[i, j])
    if n == 1
        names = "R"
    else
        names = map(i -> "R_$i", 1:n)
    end
    R = (AbsMax = (errors = e, names = names), )
    return R
end