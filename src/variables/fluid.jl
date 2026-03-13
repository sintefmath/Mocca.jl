struct AverageMolecularMass <: Jutul.ScalarVariable end

struct Concentrations <: JutulDarcy.ComponentVariables end

struct DiffusionTransmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::DiffusionTransmissibilities) = 1e-10
Jutul.minimum_value(::DiffusionTransmissibilities) = 0.0
Jutul.default_value(model, ::DiffusionTransmissibilities) = 1e-3
Jutul.associated_entity(::DiffusionTransmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::DiffusionTransmissibilities, symb)
    if haskey(data_domain, :diffusion_coefficient, Jutul.Cells())
        U = data_domain[:diffusion_coefficient]
        g = Jutul.physical_representation(data_domain)
        T = Jutul.compute_face_trans(g, model.system.p.Φ .* U)
    else
        error(":diffusion_coefficient symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

Jutul.@jutul_secondary function update_total_concentration!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G, S}, Pressure, Temperature, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    for cell in ix
        ctot[cell] = Pressure[cell] / (sys.p.R * Temperature[cell])
    end
end

Jutul.@jutul_secondary function update_average_molecular_mass!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    for cell in ix
        molecularMassOfCO2 = sys.p.molecularMassOfCO2
        molecularMassOfN2 = sys.p.molecularMassOfN2
        avm[cell] = y[CO2INDEX, cell] * molecularMassOfCO2 + y[N2INDEX, cell] * molecularMassOfN2
    end
end

Jutul.@jutul_secondary function update_concentrations!(concentrations, tv::Concentrations, model::Jutul.SimulationModel{G, S}, y, cTot, ix) where {G, S <: AdsorptionSystem}
    for cell in ix
        for component in 1:JutulDarcy.number_of_components(model.system)
            concentrations[component, cell] = y[component, cell] * cTot[cell]
        end
    end
end

Jutul.@jutul_secondary function update_total_masses!(totmass, tv::JutulDarcy.TotalMasses, model::AdsorptionModel, concentrations, FluidVolume, ix)
    for cell in ix
        for component in 1:JutulDarcy.number_of_components(model.system)
            totmass[component, cell] = concentrations[component, cell] * FluidVolume[cell]
        end
    end
end
