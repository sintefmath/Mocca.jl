
Jutul.@jutul_secondary function update_our_total_masses!(
    totmass,
    tv::JutulDarcy.TotalMasses,
    model::AdsorptionModel,
    concentrations,
    FluidVolume,
    ix
)
    sys = model.system
    for cell in ix
        totmass[1, cell] = concentrations[1, cell] * FluidVolume[cell]
        totmass[2, cell] = concentrations[2, cell] * FluidVolume[cell]
    end
end

Jutul.@jutul_secondary function update_adsorption_mass_transfer(
    adsorption_mass_transfer,
    tv::AdsorptionMassTransfer,
    model::AdsorptionModel,
    concentrations,
    Temperature,
    AdsorbedConcentration,
    ix
)
    N = JutulDarcy.number_of_components(model.system) # Statically known.
    for cell in ix
        C = @view concentrations[:, cell]
        qstar = compute_equilibrium(model.system, C, Temperature[cell])
        k = compute_ki(model.system, C, qstar)
        for i in 1:N
            adsorption_mass_transfer[i, cell] = k[i] * (qstar[i] - AdsorbedConcentration[i, cell])
        end
    end
end

Jutul.@jutul_secondary function update_column_conserved_energy(column_energy, tv::ColumnEnergy, model::Jutul.SimulationModel{G,S}, SolidVolume, C_pa, ΔH, AdsorbedConcentration, Temperature, C_pg, Pressure, avm, ix) where {G,S<:AdsorptionSystem}
    for cx in ix
        column_energy[cx] = Temperature[cx]
    end
end

Jutul.@jutul_secondary function update_wall_conserved_energy(wall_energy, tv::WallEnergy, model::Jutul.SimulationModel{G,S}, WallTemperature, ix) where {G,S<:AdsorptionSystem}
    for cellindex in ix
        wall_energy[cellindex] = WallTemperature[cellindex]
    end
end





function compute_equilibrium(sys::AdsorptionSystem, concentration, temperature)
    N = JutulDarcy.number_of_components(sys)
    T = eltype(concentration)
    qstar = @MVector zeros(T, N) # TODO: Use svector
    b = @MVector zeros(T, N) # TODO: Use svector
    d = @MVector zeros(T, N) # TODO: Use svector
    bC_sum = zero(T)
    dC_sum = zero(T)
    @inbounds for i in 1:JutulDarcy.number_of_components(sys)
        b_i = sys.p.b0[i] * exp(-sys.p.ΔUbi[i] / (sys.p.R * temperature))
        d_i = sys.p.d0[i] * exp(-sys.p.ΔUdi[i] / (sys.p.R * temperature))

        b[i] = b_i
        d[i] = d_i

        bC_sum += b_i*concentration[i]
        dC_sum += d_i*concentration[i]
    end

    @inbounds for i in 1:N
        qstar[i] = sys.p.qsbi[i] * b[i] * concentration[i] / (one(T) + bC_sum) + sys.p.qsdi[i] * d[i] * concentration[i] / (one(T) + dC_sum)
    end
    return SVector{N, T}(qstar)
end

function compute_ki(sys::AdsorptionSystem, concentration, qstar)
    D_p = sys.p.D_m / sys.p.τ
    r_p = sys.p.d_p / 2.0

    return concentration ./ qstar .* 15 .* sys.p.ϵ_p .* D_p ./ r_p^2
end

Jutul.@jutul_secondary function update_cTot!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G,S}, Pressure, Temperature, ix) where {G,S<:AdsorptionSystem}
    # Update cTot
    sys = model.system

    for cellindex in ix
        ctot[cellindex] = Pressure[cellindex] / (sys.p.R * Temperature[cellindex])
    end
end

Jutul.@jutul_secondary function update_avm!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionSystem}
    #println("Updating avm")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        sys = model.system
        molecularMassOfCO2 = sys.p.molecularMassOfCO2
        molecularMassOfN2 = sys.p.molecularMassOfN2
        avm[cellindex] = y[CO2INDEX, cellindex] * molecularMassOfCO2 + y[N2INDEX, cellindex] * molecularMassOfN2
    end
end

Jutul.@jutul_secondary function update_concentrations!(concentrations, tv::Concentrations, model::Jutul.SimulationModel{G,S}, y, cTot, ix) where {G,S<:AdsorptionSystem}
    # println("Updating concentrations")
    for cellindex in ix
        # TODO: Fix masses such that they are single parameters for the whole grid
        for component in 1:JutulDarcy.number_of_components(model.system)
            concentrations[component, cellindex] = y[component, cellindex] * cTot[cellindex]
        end
    end
end





Jutul.@jutul_secondary function update_enthalpy_change(ΔH, tv::EnthalpyChange, model::Jutul.SimulationModel{G,S}, AdsorbedConcentration, ix) where {G,S<:AdsorptionSystem}
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

Jutul.@jutul_secondary function update_cpa(cpa, tv::SpecificHeatCapacityAdsorbent, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionSystem}
    sys = model.system
    for cellindex in ix
        cpa[cellindex] = sum(y[:, cellindex] .* sys.p.C_pa)
    end
end

Jutul.@jutul_secondary function update_cpg(cpg, tv::SpecificHeatCapacityFluid, model::Jutul.SimulationModel{G,S}, y, ix) where {G,S<:AdsorptionSystem}
    sys = model.system
    for cellindex in ix
        cpg[cellindex] = sum(y[:, cellindex] .* sys.p.C_pg)
    end
end
