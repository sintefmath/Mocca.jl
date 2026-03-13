struct AdsorptionMassTransfer <: JutulDarcy.ComponentVariables end

Jutul.@jutul_secondary function update_adsorption_mass_transfer!(
    adsorption_mass_transfer,
    tv::AdsorptionMassTransfer,
    model::AdsorptionModel,
    concentrations,
    Temperature,
    AdsorbedConcentration,
    ix
)
    sys = model.system
    iso = sys.isotherm
    mt = sys.mass_transfer
    N = JutulDarcy.number_of_components(sys)
    T = eltype(adsorption_mass_transfer)
    for cell in ix
        C = SVector{N, T}(@view concentrations[:, cell])
        q = SVector{N, T}(@view AdsorbedConcentration[:, cell])
        qstar = compute_equilibrium(iso, C, Temperature[cell])
        rate = compute_mass_transfer_rate(mt, C, q, qstar)
        for i in 1:N
            adsorption_mass_transfer[i, cell] = rate[i]
        end
    end
end
