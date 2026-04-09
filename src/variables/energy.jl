abstract type Energy <: Jutul.ScalarVariable end
struct ColumnEnergy <: Energy end
struct WallEnergy <: Energy end

struct EnthalpyChange <: JutulDarcy.ComponentVariables end

abstract type SpecificHeatCapacity <: Jutul.ScalarVariable end
struct SpecificHeatCapacityAdsorbent <: SpecificHeatCapacity end
struct SpecificHeatCapacityFluid <: SpecificHeatCapacity end

struct ThermalConductivities <: Jutul.ScalarVariable end
Jutul.variable_scale(::ThermalConductivities) = 1e-10
Jutul.minimum_value(::ThermalConductivities) = 0.0
Jutul.default_value(model, ::ThermalConductivities) = 1e-3
Jutul.associated_entity(::ThermalConductivities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalConductivities, symb)
    if haskey(data_domain, :thermal_conductivity, Column())
        K = first(data_domain[:thermal_conductivity, Column()])
        g = Jutul.physical_representation(data_domain)
        nc = Jutul.number_of_cells(g)
        T = Jutul.compute_face_trans(g, fill(K, nc))
    else
        error(":thermal_conductivity on Column() must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

Jutul.@jutul_secondary function update_column_conserved_energy!(column_energy, tv::ColumnEnergy, model::Jutul.SimulationModel{G, S}, Temperature, ix) where {G, S <: AdsorptionSystem}
    for cell in ix
        column_energy[cell] = Temperature[cell]
    end
end

Jutul.@jutul_secondary function update_wall_conserved_energy!(wall_energy, tv::WallEnergy, model::Jutul.SimulationModel{G, S}, WallTemperature, ix) where {G, S <: AdsorptionSystem}
    for cell in ix
        wall_energy[cell] = WallTemperature[cell]
    end
end

Jutul.@jutul_secondary function update_enthalpy_change!(ΔH, tv::EnthalpyChange, model::Jutul.SimulationModel{G, S}, concentrations, Temperature, ix) where {G, S <: AdsorptionSystem}
    iso = model.system.isotherm
    N = JutulDarcy.number_of_components(model.system)
    T = eltype(ΔH)
    for cell in ix
        C = SVector{N, T}(@view concentrations[:, cell])
        ΔH_values = compute_enthalpy(iso, C, Temperature[cell])
        for i in 1:N
            ΔH[i, cell] = ΔH_values[i]
        end
    end
end

Jutul.@jutul_secondary function update_heat_capacity_adsorbent!(cpa, tv::SpecificHeatCapacityAdsorbent, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    C_pa = sys.heat_capacity_adsorbed
    T = eltype(cpa)
    for cell in ix
        cpa_i = zero(T)
        for component in 1:JutulDarcy.number_of_components(sys)
            cpa_i += y[component, cell] * C_pa[component]
        end
        cpa[cell] = cpa_i
    end
end

Jutul.@jutul_secondary function update_heat_capacity_fluid!(cpg, tv::SpecificHeatCapacityFluid, model::Jutul.SimulationModel{G, S}, y, ix) where {G, S <: AdsorptionSystem}
    sys = model.system
    C_pg = sys.heat_capacity_gas
    T = eltype(cpg)
    for cell in ix
        cpg_i = zero(T)
        for component in 1:JutulDarcy.number_of_components(sys)
            cpg_i += y[component, cell] * C_pg[component]
        end
        cpg[cell] = cpg_i
    end
end
