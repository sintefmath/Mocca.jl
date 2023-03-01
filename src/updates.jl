
Jutul.@jutul_secondary function update_our_total_masses!(
    totmass,
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    concentrations,
    ix
)
    sys = model.system
    pv = Jutul.physical_representation(model.domain).pore_volumes
    # pv = sys.Φ * ones(size(pv)) # FIXME: Remove this
    for cell in ix
        #@info "Data in cell" cTot[cell] y[:, cell] PhaseMassDensities[1, cell]  model.domain.grid.pore_volumes[cell]
        totmass[1, cell] = concentrations[1, cell] * pv[cell]
        totmass[2, cell] = concentrations[2, cell] * pv[cell]
        # totmass[CO2COMPONENTINDEX, cell] = y
        # totmass[N2COMPONENTINDEX, cell] = y
    end
end

Jutul.@jutul_secondary function update_adsorption_mass_transfer(
    adsorption_mass_transfer,
    tv::AdsorptionMassTransfer,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    concentrations,
    Temperature,
    adsorptionRates,
    ix
)
    for cell in ix
        qstar = compute_equilibrium(model.system, concentrations[:, cell], Temperature[cell])
        k = compute_ki(model.system, concentrations[:, cell], qstar)
        force = k .* (qstar .- adsorptionRates[:, cell])

        # @info "cell $ix" qstar k concentrations adsorptionRates force

        adsorption_mass_transfer[:, cell] = force
    end
end

function compute_equilibrium(sys::AdsorptionFlowSystem, concentration, temperature)
    qstar = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    b = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    d = zeros(eltype(concentration), JutulDarcy.number_of_components(sys)) # TODO: Use svector
    for i in 1:JutulDarcy.number_of_components(sys)
        b[i] = sys.b0[i] * exp(-sys.ΔUbi[i] / (sys.R * temperature))
        d[i] = sys.d0[i] * exp(-sys.ΔUdi[i] / (sys.R * temperature))
    end

    for i in 1:JutulDarcy.number_of_components(sys)
        qstar[i] = sys.qsbi[i] * b[i] * concentration[i] / (1 + sum(b .* concentration)) + sys.qsdi[i] * d[i] * concentration[i] / (1 + sum(d .* concentration))
    end
    return qstar
end

function compute_ki(sys::AdsorptionFlowSystem, concentration, qstar)
    D_p = sys.D_m / sys.τ
    r_p = sys.d_p / 2.0

    return concentration ./ qstar .* 15 .* sys.ϵ_p .* D_p ./ r_p^2
end

Jutul.@jutul_secondary function update_cTot!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G,S}, Pressure, Temperature, ix) where {G,S<:AdsorptionFlowSystem}
    # Update cTot
    sys = model.system

    for cellindex in ix
        ctot[cellindex] = Pressure[cellindex] / (sys.R * Temperature[cellindex])
    end
end

Jutul.@jutul_secondary function update_avm!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionFlowSystem}
    #println("Updating avm")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        sys = model.system
        molecularMassOfCO2 = sys.molecularMassOfCO2
        molecularMassOfN2 = sys.molecularMassOfN2
        avm[cellindex] = y[CO2INDEX, cellindex] * molecularMassOfCO2 + y[N2INDEX, cellindex] * molecularMassOfN2
    end
end

Jutul.@jutul_secondary function update_concentrations!(concentrations, tv::Concentrations, model::Jutul.SimulationModel{G,S}, y, cTot, ix) where {G,S<:AdsorptionFlowSystem}
    # println("Updating concentrations")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        for component in 1:JutulDarcy.number_of_components(model.system)
            concentrations[component, cellindex] = y[component, cellindex] * cTot[cellindex]
        end
    end
end
