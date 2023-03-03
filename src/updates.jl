
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
        b[i] = sys.p.b0[i] * exp(-sys.p.ΔUbi[i] / (sys.p.R * temperature))
        d[i] = sys.p.d0[i] * exp(-sys.p.ΔUdi[i] / (sys.p.R * temperature))
    end

    for i in 1:JutulDarcy.number_of_components(sys)
        qstar[i] = sys.p.qsbi[i] * b[i] * concentration[i] / (1 + sum(b .* concentration)) + sys.p.qsdi[i] * d[i] * concentration[i] / (1 + sum(d .* concentration))
    end
    return qstar
end

function compute_ki(sys::AdsorptionFlowSystem, concentration, qstar)
    D_p = sys.p.D_m / sys.p.τ
    r_p = sys.p.d_p / 2.0

    return concentration ./ qstar .* 15 .* sys.p.ϵ_p .* D_p ./ r_p^2
end

Jutul.@jutul_secondary function update_cTot!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G,S}, Pressure, Temperature, ix) where {G,S<:AdsorptionFlowSystem}
    # Update cTot
    sys = model.system

    for cellindex in ix
        ctot[cellindex] = Pressure[cellindex] / (sys.p.R * Temperature[cellindex])
    end
end

Jutul.@jutul_secondary function update_avm!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionFlowSystem}
    #println("Updating avm")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        sys = model.system
        molecularMassOfCO2 = sys.p.molecularMassOfCO2
        molecularMassOfN2 = sys.p.molecularMassOfN2
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



Jutul.@jutul_secondary function update_column_conserved_energy(column_energy, tv::ColumnEnergy, model::Jutul.SimulationModel{G,S}, solidVolume, C_pa, avm, adsorptionRates, Temperature, ix) where {G,S<:AdsorptionFlowSystem}
    sys = model.system
    ρ_s = sys.p.ρ_s
    C_ps = sys.p.C_ps
    for cellindex in ix
        sq = sum(adsorptionRates[:, cellindex])
        column_energy[cellindex] = solidVolume[cellindex] * (ρ_s * C_ps + C_pa[cellindex] * avm[cellindex] * sq) * Temperature[cellindex]
    end
end

Jutul.@jutul_secondary function update_wall_conserved_energy(wall_energy, tv::WallEnergy, model::Jutul.SimulationModel{G,S}, WallTemperature, ix) where {G,S<:AdsorptionFlowSystem}
    sys = model.system
    for cellindex in ix
        wall_energy[cellindex] = sys.p.ρ_w * sys.p.C_pw * WallTemperature[cellindex]
    end
end

Jutul.@jutul_secondary function update_enthalpy_change(ΔH, tv::EnthalpyChange, model::Jutul.SimulationModel{G,S}, adsorptionRates, ix) where {G,S<:AdsorptionFlowSystem}
    sys = model.system

    qsbi = sys.p.qsbi
    qsdi = sys.p.qsdi
    sumq = qsbi[CO2INDEX] + qsdi[CO2INDEX]
    ΔUbi = sys.p.ΔUbi
    ΔUdi = sys.p.ΔUdi
    R = sys.p.R
    T0 = sys.p.T0 # TODO: Review
    for cellindex in ix
        for i in eachindex(ΔH[:, cellindex])
            ΔH[i, cellindex] = (qsbi[i] * (ΔUbi[i] - R * T0) + qsdi[i] * (ΔUdi[i] - R * T0)) / sumq
        end
    end
end

Jutul.@jutul_secondary function update_cpa(cpa, tv::SpecificHeatCapasityAdsorbent, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionFlowSystem}
    sys = model.system
    for cellindex in ix
        cpa[cellindex] = sum(y[:, cellindex] .* sys.p.C_pa)
    end
end

Jutul.@jutul_secondary function update_cpg(cpg, tv::SpecificHeatCapasityFluid, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionFlowSystem}
    sys = model.system
    for cellindex in ix
        cpg[cellindex] = sum(y[:, cellindex] .* sys.p.C_pg)
    end
end
