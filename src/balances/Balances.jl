"""
    Mocca.Balances

The material and energy balances shared by every unit: gas, adsorbed loading,
gas energy and wall energy, with the source terms that enter them, their
variables and parameters, scale-free convergence criteria, and the setup
shared by every bed of sorbent.
"""
module Balances

import Jutul
using StaticArrays
using ..Framework, ..Physics
using ..Framework: _gas_state, _state_has, _state_or_domain

export AdsorbedConcentration, AdsorbedLoading, AdsorbedMolarHeatCapacity,
    AdsorbedPhaseHeatCapacity, AdsorbentDensity, AdsorbentHeatCapacity,
    AdsorptionMassTransfer, AmbientLoss, AmbientTemperature, AverageMolarMass,
    AxialDispersion, BedCrossSectionArea, BulkVolume, calc_bc_trans, calc_bc_wall_trans,
    CellDx, column_mesh, ColumnEnergy, ComponentMasses, compute_column_face_area,
    ConservedPressure, DispersionScale, DispersionTransmissibilities, EnthalpyChange,
    FluidViscosity, FluidVolume, GasMolarHeatCapacity, HeatOfSorption,
    InnerHeatTransferCoeff, mocca_domain, MolarConcentration, MolecularDiffusivity,
    OuterHeatTransferCoeff, Permeability, ScaledAdsorptionMassTransfer,
    setup_process_model, setup_process_parameters, setup_process_state, SolidHeatCapacity,
    SorptionUptake, ThermalConductivities, TimeScale, TotalMolarConcentration,
    ViscosityFactor, WallArea, WallConductivity, WallCrossSectionArea, WallDensity,
    WallEndConduction, WallEnergy, WallExchange, WallHeatCapacity

include("geometry.jl")
include("unit_parameters.jl")
include("gas.jl")
include("sorbent.jl")
include("energy.jl")
include("wall.jl")
include("convergence.jl")
include("sorbent_bed.jl")

end
