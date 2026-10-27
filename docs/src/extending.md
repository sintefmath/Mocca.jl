```@meta
CurrentModule = Mocca
```

# Extending Mocca

Most studies need only the [Reference](@ref) functions. This page is for adding
your own kinds of step and unit, and for working with the simulation model
that a [`Process`](@ref) builds. That model is a
[Jutul](https://github.com/sintefmath/Jutul.jl) `MultiModel`, with one model
per bed and device, and the steps of the cycle become its forces.

## Adding a kind of process

A kind of process is a subtype of [`Process`](@ref), and its result one of
[`ProcessResult`](@ref) with fields `process`, `states` and `timesteps`. It
needs fields `flowsheet`, `options` ([`ProcessOptions`](@ref)), `model`,
`parameters`, `parameter_settings`, `variables`, `integrals` and `measures`,
and these methods:

- [`simulate`](@ref), [`initial_state`](@ref) and [`steady_state`](@ref),
  which finds the state it settles to;
- [`stages`](@ref), the device settings and durations it is operated by;
- [`performance_quantities`](@ref) and [`performance`](@ref) of its result,
  from [`performance_sums`](@ref) over the steps that make up one period of
  its operation;
- for design variables of its own kinds of place, `_get(p, place)` and
  `_set!(p, place, value)` (a [`ParameterPlace`](@ref) works for every kind),
  and `_rebuild!(p)`, which sets up its model and parameters again after a
  [`GeometryPlace`](@ref) rebuilt a unit; gradients then need nothing more;
- for gradients, `_fixed_point(p, result)`: the forces and time steps of one
  period, the state at its start, at the fixed point the result settled to,
  and the result to measure.

[`performance_integrals`](@ref), [`performance_measures`](@ref),
[`add_measure!`](@ref), [`add_variable!`](@ref), [`performance`](@ref) of
the process, [`performance_gradient`](@ref) and [`optimise!`](@ref) then
work unchanged.

```@docs
ProcessOptions
```

## Adding a kind of place

A new kind of place, such as a field of a device or unit that no place
reaches yet, is a subtype of [`Place`](@ref) with methods `_get(p, place)` and
`_set!(p, place, value)` for the kinds of process that have it, and
`_describe(place)` for messages. Its gradients then come from the general
design Jacobian, with no derivative code. To read and write it in process
files, add methods of `_place_to_dict` and a case in `_place_from_dict`.

## Adding a kind of step

The set of steps is open: define a subtype of [`AbstractStep`](@ref) and a
method of [`step_settings`](@ref).

```@docs
AbstractStep
StepContext
compile_cycle
step_settings
end_pressure
step_beds
default_step_name
register_step!
find_device
valve
boundary
start_pressure
ramp
```

## Adding a thermodynamics model

The heat capacities and the equation of state of the gas come from an
[`AbstractGasThermodynamics`](@ref) model, held by the [`Gas`](@ref) and by
every unit with a gas phase (the `thermodynamics` field of its system; a bed
also has `adsorbed_thermodynamics` for the adsorbed phase). Mocca has
[`IdealGasConstantCp`](@ref) and [`IdealGasShomate`](@ref). Another model,
such as NASA polynomials, is a subtype with one method,
[`component_heat_capacities`](@ref): the molar heat capacity of each
component [J/(mol·K)] at a temperature, as an `SVector`.

```julia
using Mocca, StaticArrays

# NASA 7-coefficient polynomials, one temperature range: cp/R = a1 + a2 T + ... + a5 T⁴
struct NASA7{N} <: Mocca.AbstractGasThermodynamics
    a::NTuple{5, SVector{N, Float64}}
end

function Mocca.component_heat_capacities(m::NASA7, T)
    a = m.a
    return Mocca.GAS_CONSTANT .* (a[1] .+ T .* (a[2] .+ T .* (a[3] .+ T .* (a[4] .+ T .* a[5]))))
end

gas = Gas(["CO2", "N2"]; thermodynamics = NASA7((SVector(2.36, 3.53), SVector(8.98e-3, -1.24e-4),
    SVector(-7.12e-6, -5.03e-7), SVector(2.46e-9, 2.44e-9), SVector(-1.44e-13, -1.41e-12))))
```

The rest follows:
[`component_enthalpies`](@ref), relative to 298.15 K, by Gauss–Legendre
quadrature of the heat capacities (add a method if the integral is known),
and from these [`molar_heat_capacity`](@ref), [`molar_enthalpy`](@ref) and
[`heat_capacity_ratio`](@ref) of a mixture. [`molar_concentration`](@ref) is
that of the ideal gas; the energy balance of the column also assumes an ideal
gas, so a model with another equation of state needs more than a method for
it.

The functions are evaluated in every cell at every Newton iteration with
automatic-differentiation numbers for `T`, so they must accept any `Real`, and
should not allocate: keep the parameters in `SVector`s or tuples, and the
struct `isbits`, as [`DualSiteLangmuir`](@ref) does.

The model's heat capacities enter the energy balance through the gas and
adsorbed-phase heat capacity variables `C_pg` and `C_pa`, the sensible heat of
sorption and the enthalpy a stream brings into a bed (see [Model
equations](@ref)). A model that is not one of Mocca's cannot be written to a
process file.

## Adding a transport model

How the viscosity and molecular diffusivities change through a bed comes from
an [`AbstractGasTransport`](@ref) model, held by the [`Gas`](@ref) and, with
its values at the feed, by the bed (a [`TransportReference`](@ref) in the
`transport` field of its system). Mocca has [`ConstantTransport`](@ref), the
default, and [`ChapmanEnskogTransport`](@ref). Another model gives two
functions: [`mixture_viscosity`](@ref) and [`molecular_diffusivities`](@ref),
the diffusivity of each component in the mixture. The bed scales the values at
the feed by their ratios, so a model may be exact in its trends and looser in
its absolute values. As for thermodynamics, they must take
automatic-differentiation numbers and not allocate.

```julia
using Mocca, StaticArrays

# Viscosity as a power of temperature, the same for every mixture, and
# Knudsen diffusion, as the square root of temperature whatever the pressure
struct PowerLawTransport <: Mocca.AbstractGasTransport
    a::Float64   # exponent of the viscosity
end
Mocca.mixture_viscosity(m::PowerLawTransport, y, T) = 1.8e-5 * (T / 298.15)^m.a
Mocca.molecular_diffusivities(m::PowerLawTransport, y, T, P) = 1.6e-5 * sqrt(T / 298.15) .* SVector(1.0, 1.0)

gas = Gas(["CO2", "N2"]; transport = PowerLawTransport(0.7))
```

A bed with a transport model other than `ConstantTransport` adds a per-cell
variable `ViscosityFactor`, and evaluates its mass transfer as
`ScaledAdsorptionMassTransfer`, which passes the diffusivity factors at the
feed pressure to [`compute_mass_transfer_rate`](@ref) after the temperature.
A mass transfer model whose rate depends on the molecular diffusivity adds a
method with that sixth argument, as [`LinearDrivingForce`](@ref) does.
Unless a `dispersion` is given,
its axial dispersion is [`WakaoFunazkri`](@ref), with a per-cell
`MolecularDiffusivity` and a `DispersionScale` parameter. A model of your own
cannot be written to a process file.

## The layer underneath

The objects of the [Reference](@ref) build these, which can also be used
directly. A
[`Flowsheet`](@ref) holds a [`UnitNetwork`](@ref): the unit models and the
connections between their ports, from which the simulation model is built.

```@docs
fixed_bed_model
MoccaSystem
DistributedUnit
LumpedUnit
FlowElement
component_names
AbstractFlowsheet
UnitNetwork
network
add_unit!
connect!
set_boundary!
PortRef
AbstractConnection
connection_ports
check_connection
add_cross_terms!
StreamConnection
holdup_units
```

## Adding a kind of unit

A new kind of unit is a `MoccaSystem` with these methods.

```@docs
stream_equations
setup_unit_state
setup_unit_parameters
holds_sorbent
AxialDispersion
```

## Rules for differentiable models

Every gradient in Mocca comes from Jutul's automatic differentiation and
adjoints of the model's own equations, so a new kind of unit, device law,
objective or measure needs no derivative code. It does need to follow these
rules, each learned the hard way; a model that breaks one gives gradients that
are wrong without any error.

1. **A branch must not change which cells, entities or parameters a term
   reads.** Jutul finds what each equation, cross term and objective depends
   on by tracing it once, at the state it starts from, often a uniform one,
   and keeps only those dependencies for the whole run. A term that reads a
   neighbour only where a slope is steep, or one port only when the gas flows
   one way, is traced along one branch and its other dependencies are
   dropped. Read every value on every branch, with zero weight where it does
   not count: `ifelse`, which evaluates both sides, rather than `if`; the
   limiters read both neighbours through `_stencil_anchor`
   (`src/balances/limiters.jl`), and cross terms touch their target through
   `_target_anchor` (`src/networks/ports.jl`). Branching between variables of
   one entity is fine, such as the upstream port of a device, whose variables
   all live in its one cell.
2. **Do not compare a value that carries derivatives with zero.** A dual
   number of value 0 with a positive derivative compares as greater than 0,
   and in an adjoint the parameters are dual numbers, so a test such as
   `Σc > 0` takes the wrong branch exactly where the amount vanishes, and
   divides by zero. Use smooth forms, or a positive threshold well away from
   the values that occur.
3. **Keep denominators, and their derivatives, finite without a branch.**
   Adding a tiny amount works for a denominator that is a sum of amounts of
   ordinary size, such as `Σc + 1e-20` for the gas in the pores
   (`src/balances/solid_advection.jl`). Where the denominator is the square of
   a small difference, as in a slope limiter, an amount small enough not to
   matter has a square that underflows in the derivative: return zero below a
   cutoff on the value instead (`_LIMITER_TINY = 1e-75`).
4. **What a design changes must reach the equations as a parameter.** A
   quantity computed once when the model is built, from the geometry or a
   flow, keeps its value when the design changes, and its effect is lost from
   the gradient. Cell sizes, volumes, areas and transmissibilities are
   parameters, and the default valve conductances are recomputed when the
   cycle compiles, so geometry gets gradients (see [`GeometryPlace`](@ref));
   a heat transfer coefficient computed from the sorbent velocity outside the
   unit does not follow it.
5. **Objectives and measures are generic Julia functions.** An integral (see
   [`add_integral!`](@ref)) reads the state and the unit's parameters, which
   carry dual numbers in the adjoints, and a measure (see
   [`add_measure!`](@ref)) is differentiated by ForwardDiff, so neither may
   convert its arguments to `Float64`.

[`check_gradients`](@ref) compares the gradients of a process with central
differences at two step sizes and reports a failure when the two differences
agree with each other but not with the gradient. `test/gradient_checks.jl`
runs it on every process of the library and every process the examples
build; set `MOCCA_GRADIENT_CHECKS=true` to include it in the tests.

```@docs
check_gradients
```

## Moving beds

A [`MovingBed`](@ref) carries its sorbent down through the unit.
[`add_moving_bed!`](@ref) builds one with [`moving_bed_model`](@ref). Sections
of a moving bed process pass sorbent to each other through a
[`SorbentConnection`](@ref) and, if they have a [`HeatExchanger`](@ref), heat
transfer fluid through a [`HeatTransferFluidConnection`](@ref). The
[Moving bed TSA](examples/mbtsa_mondino_2022.md) example builds the process of
Mondino et al. (2022). What the sorbent and the fluid carry is reconstructed at
the cell faces with a slope limiter; a new limiter is a subtype of
[`AbstractLimiter`](@ref) with a method of [`limited_slope`](@ref).

```@docs
MovingBed
setup_moving_bed_model
tube_bundle
moving_bed_sorbent_htc
HeatExchanger
moves_sorbent
SorbentConnection
connect_solids!
HeatTransferFluidConnection
connect_heat_transfer_fluid!
sorbent_path
AbstractLimiter
limited_slope
Upwind
VanAlbada
VanLeer
Minmod
advection_limiter
SolidStreamCT
HXFluidStreamCT
FixedMolarFlow
BackPressureController
Virial
moving_bed_model
MovingBedSpec
stream_rates
heat_duty
equivalent_work
hx_fluid_inlet_temperature
moving_sorbent_units
sorbent_inventory
sorbent_circulation_time
plot_moving_bed_profiles
```

## Devices, schedules and steady states

```@docs
FlowDevice
setup_device_model
AbstractDeviceLaw
device_flow
LinearValve
VolumetricFlow
PortStateCT
StreamCT
Stage
setup_schedule
parallel_stages
PortCondition
ExponentialRamp
port_conductance
default_conductance
setup_flowsheet_model
setup_flowsheet_state
setup_flowsheet_parameters
DualSiteLangmuir
AbstractIsotherm
AbstractMassTransfer
compute_mass_transfer_rate
LinearDrivingForce
pore_diffusivity
has_knudsen
simulate_to_cyclic_steady_state
newton_cyclic_steady_state
simulate_until_mass_balance
CycleProgress
cycle_change
cyclic_steady_state_gradient
force_gradients
StreamObjective
VacuumPumpObjective
CompressorObjective
BlowerObjective
HeatDutyObjective
EquivalentWorkObjective
SumObjective
objective_state
integral_values
measure_partials
SimulationFailure
```
