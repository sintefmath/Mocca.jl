"""
    Mocca.Balances

The material and energy balances shared by every unit: gas, adsorbed loading,
gas energy and wall energy, sorbent carried through a unit, and tubes of heat
transfer fluid, with the source terms that enter them, their
variables and parameters, scale-free convergence criteria, and the setup
shared by every bed of sorbent.
"""
module Balances

import Jutul
using StaticArrays
using ..Framework, ..Physics
using ..Framework: _gas_state, _state_has, _state_or_domain

export AbstractLimiter, advection_limiter, limited_slope, Minmod, Upwind, VanAlbada, VanLeer
export AdsorbedConcentration, AdsorbedLoading, AdsorbedMolarHeatCapacity,
    AdsorbedPhaseHeatCapacity, AdsorbentDensity, AdsorbentHeatCapacity,
    AdsorptionMassTransfer, AmbientLoss, AmbientTemperature, AverageMolarMass,
    AxialDispersion, BedCrossSectionArea, BulkVolume, calc_bc_trans, calc_bc_wall_trans,
    CellDx, column_mesh, ColumnEnergy, ComponentMasses, compute_column_face_area,
    compute_dispersion, compute_permeability, ConservedPressure, DispersionScale,
    DispersionTransmissibilities, EnthalpyChange, FluidSideHeatTransferCoeff, FluidViscosity,
    FluidVolume, GasMolarHeatCapacity, HeatOfSorption, hx_fluid_heat_capacity_flow,
    HXFluidDensity, HXFluidEnergy, HXFluidHeatCapacity, HXFluidInletTemperature,
    HXFluidMassFlow, HXFluidVolume, InnerHeatTransferCoeff, mocca_domain, MolarConcentration,
    MolecularDiffusivity, OuterHeatTransferCoeff, ParticlePorosity, Permeability,
    pore_gas_energy, pore_gas_volume_flow, PoreGasFlow, ScaledAdsorptionMassTransfer,
    select_heat_exchanger_parameters!, setup_process_model, setup_process_parameters,
    setup_process_state, SolidHeatCapacity, SolidInletGasConcentration, SolidInletLoading,
    SolidInletTemperature, sorbent_heat_capacity_flow, sorbent_volume_flow, SorbentHeatFlow,
    SorbentMassFlow, SorbentSideHeatTransferCoeff, SorptionRate, SorptionRateHeat,
    SorptionUptake, ThermalConductivities, TimeScale, TotalMolarConcentration, TubeDensity,
    TubeEnergy, TubeExchange, TubeHeatCapacity, TubeInnerArea, TubeOuterArea, TubeVolume,
    ViscosityFactor, WallArea, WallConductivity, WallCrossSectionArea, WallDensity,
    WallEndConduction, WallEnergy, WallExchange, WallHeatCapacity

include("geometry.jl")
include("unit_parameters.jl")
include("gas.jl")
include("sorbent.jl")
include("energy.jl")
include("wall.jl")
include("limiters.jl")
include("solid_advection.jl")
include("heat_exchanger.jl")
include("convergence.jl")
include("correlations.jl")
include("sorbent_bed.jl")

end
