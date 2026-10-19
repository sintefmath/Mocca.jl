"""
    Mocca.Framework

System hierarchy, entities, core variables, and the interface every unit
implements.
"""
module Framework

import Jutul
import LinearAlgebra

export AbstractSourceTerm, AbstractThermalModel, Adiabatic, component_names,
    ComponentVariable, DistributedUnit, DistributedUnitModel, energy_source,
    energy_source_terms, FixedTemperatureEquation, FlowChannel, FlowElement, GAS_CONSTANT,
    GasMoleFractions, has_energy_balance, has_wall, heat_capacity, heat_capacity_terms,
    holds_sorbent, Isothermal, IsothermalTemperature, LumpedUnit, mass_source,
    mass_source_terms, MoccaCase, MoccaModel, MoccaSystem, number_of_components, ports,
    Pressure, setup_unit_parameters, setup_unit_state, SimulationFailure, SorbentBed,
    SorbentBedModel, stream_equations, sum_terms, Temperature, thermal_model, Unit,
    wall_source, wall_source_terms, WithWall

include("types.jl")
include("variables.jl")
include("conservation.jl")
include("thermal.jl")
include("unit_interface.jl")

end
