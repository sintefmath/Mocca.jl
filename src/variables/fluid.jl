struct AverageMolecularMass <: Jutul.ScalarVariable end

struct Concentrations <: JutulDarcy.ComponentVariables end

struct DiffusionTransmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::DiffusionTransmissibilities) = 1e-10
Jutul.minimum_value(::DiffusionTransmissibilities) = 0.0
Jutul.default_value(model, ::DiffusionTransmissibilities) = 1e-3
Jutul.associated_entity(::DiffusionTransmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::DiffusionTransmissibilities, symb)
    if haskey(data_domain, :diffusion_coefficient, Column())
        D = first(data_domain[:diffusion_coefficient, Column()])
        g = Jutul.physical_representation(data_domain)
        porosity = data_domain[:porosity]
        T = Jutul.compute_face_trans(g, porosity .* D)
    else
        error(":diffusion_coefficient on Column() must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

Jutul.@jutul_secondary function update_total_concentration!(ctot, tv::JutulDarcy.TotalMass, model::Jutul.SimulationModel{G, S}, Pressure, Temperature, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    for cell in ix
        ctot[cell] = Pressure[cell] / (GAS_CONSTANT * Temperature[cell])
    end
end

Jutul.@jutul_secondary function update_average_molecular_mass!(avm, tv::AverageMolecularMass, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    mm = sys.molecular_masses
    for cell in ix
        avm_val = zero(eltype(avm))
        for component in 1:JutulDarcy.number_of_components(sys)
            avm_val += y[component, cell] * mm[component]
        end
        avm[cell] = avm_val
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
