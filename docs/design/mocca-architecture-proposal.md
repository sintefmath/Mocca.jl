# Mocca as a modular, differentiable process simulator

*Architecture proposal — draft for discussion. 25 September 2026. Based on Mocca.jl `main` at v0.1.0 (commit `70f6cd6`) and Jutul 0.4.31.*

*Status, 1 October 2026: phases 0–2 and 4 are implemented, 3 in part. Section 10 lists what was built, where it departs from this proposal, and issues found in the existing code.*

Mocca simulates one fixed-bed adsorption column today. This proposal restructures it into a set of unit models — fixed, moving, rotating and fluidised beds, membranes, solvent absorbers, and equipment such as valves, pumps and compressors — that connect into a flowsheet through Jutul `MultiModel`. The whole flowsheet stays differentiable.

## What changes

- **One set of equations for every unit.** Each unit supplies the source terms and fluxes that are specific to it.
- **Units connect through ports.** Boundary conditions that currently describe the feed or vacuum pump become separate device models, and a cycle becomes a schedule of valve openings and set-points.
- **All tunable quantities become Jutul parameters or forces**, so the adjoint can compute gradients for them.
- **Existing work is kept.** The isotherm and mass-transfer interfaces, the membrane `MultiModel` prototype and the regression tests carry over, and single-column simulations keep working.

## 1. Starting point

The proposal builds on three pieces of existing work: the fixed-bed code on `main`, the `MultiModel` membrane prototype on `claude/membrane-unit-modeling-0178f3`, and the mass and energy balance functions on `copilot/setup-mass-energy-balance-functions`.

### Worth keeping as it is

- **Swappable physics objects.** `AbstractIsotherm` and `AbstractMassTransfer` are small interfaces held as type parameters on the system, so any implementation can be plugged in and still works with automatic differentiation (AD).
- **Structs generic in `RealT`.** This is what lets the current `DictOptimization` example differentiate through case setup.
- **The membrane branch's base types.** It adds a `MoccaSystem`/`MoccaModel` supertype, `AbstractPermeanceModel` and `AbstractReactionModel`, and the first `MultiModel` in the codebase. That model couples a tube and a shell side through `AdditiveCrossTerm`s with `CTSkewSymmetry`, which applies one computed flux with opposite signs to both sides.

### What blocks a multi-unit simulator

1. **The equations are specific to one reactor.** The `:ComponentMasses` residual in `src/equations/flow/flow.jl` has the `V·∂q/∂t` sorption term written into it, so each new reactor would need its own copy of the equation.
2. **Boundary conditions also describe the neighbouring equipment.** `PressurisationBC`, `EvacuationBC` and the others are tied to cell 1 and cell N and include the feed or vacuum pump's behaviour. Two units can't be connected this way.
3. **Some values bypass the parameter system, so the adjoint can't compute gradients for them.** The boundary conditions read `fluid_viscosity`, `permeability`, `r_in` and `ambient_temperature` directly from `model.data_domain`. `FluidDensity` is also a constant, although the gas is modelled as ideal.
4. **The `Column()` entity name assumes a column.**

## 2. Layers

Code depends only on the layers below it. `L1`, `L3`, `L4` and `L5` are new; `L0` and `L2` reorganise existing code.

| Layer | Directory | Contents |
|---|---|---|
| L5 | `process/` | Cycle schedules, cyclic steady state (CSS), performance metrics, objectives, optimisation |
| L4 | `flowsheet/` | Build a `Jutul.MultiModel` from a flowsheet, `connect!`, per-unit state and parameter setup |
| L3 | `coupling/` | Ports, cross terms for streams and interfaces, counter-current pairing |
| L2 | `units/` | Fixed, moving, rotating and fluidised beds, membrane, absorber, equipment |
| L1 | `blocks/` | Reusable sets of variables, fluxes and source terms: gas flow, sorbent, wall energy, solid advection, liquid flow |
| L0 | `physics/`, `thermo/` | Pure functions dispatched on physics objects: isotherms, LDF, permeance, kinetics, equation of state (EOS), Ergun, dispersion and heat-transfer correlations |

Proposed `src/` layout: `core/` (system hierarchy, entities, traits, conservation-law skeleton, convergence), `thermo/`, `physics/{isotherms, mass_transfer, permeance, reactions, gas_liquid, transport}`, `blocks/`, `units/<type>/`, `coupling/`, `flowsheet/`, `process/`, `io/`. Plotting moves into a Makie package extension, so `CairoMakie` no longer loads with every `using Mocca`.

## 3. System type hierarchy

Every unit is a Jutul `SimulationModel` whose system sits in one abstract hierarchy:

```julia
abstract type MoccaSystem     <: Jutul.JutulSystem end   # from the membrane branch
abstract type DistributedUnit <: MoccaSystem end          # has a mesh and solves PDEs
abstract type SorbentBed      <: DistributedUnit end      # gas + sorbent in the same cell
abstract type FlowChannel     <: DistributedUnit end      # one side of a two-channel device
abstract type LumpedUnit      <: MoccaSystem end          # 0D with holdup: tank, node, flash drum
abstract type FlowElement     <: MoccaSystem end          # no holdup: valve, pump, compressor, vacuum pump
abstract type BoundaryUnit    <: MoccaSystem end          # feed, product, ambient (see section 10: implemented as port conditions)

struct FixedBed{N,R,Iso,MT,Th,Mom,Eos} <: SorbentBed  end   # today's AdsorptionSystem
struct MovingBed{...}                  <: SorbentBed  end
struct RotatingBed{...}                <: SorbentBed  end
struct FluidisedBed{...}               <: SorbentBed  end   # bubble and emulsion regions
struct GasChannel{Side,...}            <: FlowChannel end   # membrane tube/shell, absorber gas side
struct LiquidChannel{...}              <: FlowChannel end   # absorber or stripper solvent side
```

Physics choices are made by passing objects into the system, not by creating subclasses:

| Choice | Options |
|---|---|
| Isotherm, mass transfer | Existing types, plus GAB from `fran/isotherm_updates` |
| `thermal` | `Isothermal`, `Adiabatic`, `WithWall` |
| `momentum` | `Darcy`, `Ergun` |
| `eos` | `IdealGas` now; Peng–Robinson later for compression and dense-phase CO₂ |
| `reaction`, `permeance`, `film_model` | Existing reaction and permeance types; a gas–liquid film model for the absorber |

Small trait functions such as `has_energy(sys)`, `has_wall(sys)`, `solid_is_mobile(sys)` and `has_liquid(sys)` drive the `select_*!` functions. The traits decide which variables and equations exist. The update code for each variable is written once, on the broadest model type it applies to, usually `MoccaModel`.

Jutul's `CompositeSystem` was considered and isn't recommended. Gas and sorbent are tightly coupled inside each cell, and splitting them into namespaced sub-systems adds complexity without benefit.

## 4. Sharing functions

This is the core refactor. Every distributed unit uses one conservation-law residual, and each unit type supplies a tuple of source terms and fluxes:

```julia
# generic, written once for every DistributedUnit:
R_i = ∂M_i/∂t + ∇·F_i − Σ_s source(s, sys, state, state0, cell, i, Δt)

mass_sources(::FixedBed)    = (SorptionUptake(),)
mass_sources(::GasChannel)  = (ReactionSource(),)     # NoReaction drops out at compile time
energy_sources(::FixedBed)  = (HeatOfSorption(), PressureWork(), WallExchange())
flux_model(::MovingBed)     = (GasAdvectionDispersion(), SolidAdvection())
```

The membrane branch's `reaction_mass_source` already does this in a small way; this generalises it. The type of the source-term tuple is known at compile time, so the loop over it is type-stable and costs nothing at runtime.

### Where each function is defined

| Level | Dispatches on | Examples |
|---|---|---|
| Physics functions | Physics object | `compute_equilibrium`, `compute_mass_transfer_rate`, `compute_permeate_flux`, EOS, Ergun |
| Secondary variables | `MoccaModel` or trait | `TotalMolarConcentration`, `MolarConcentration`, `AverageMolarMass`, `C_pg`, gas density |
| Blocks | Source or flux type | `SorptionUptake`, `WallExchange`, `SolidAdvection`, `GasAdvectionDispersion` |
| Unit-specific | Concrete system | Rotating-bed sector geometry, fluidised-bed bubble correlations |
| Cross terms | Cross-term type | `StreamCT`, `PermeationCT`, `FilmTransferCT`, `HeatExchangeCT`, `SolidsTransferCT` |

### Which blocks each reactor reuses

| Block | Fixed | Moving | Rotating | Fluidised | Membrane | Absorber | Equipment |
|---|---|---|---|---|---|---|---|
| 1D gas flow (Darcy/Ergun + dispersion) | ● | ● | ● axial | ● per region | ● each side | ● gas side | – |
| Sorbent loading (isotherm + LDF) | ● | ● | ● | ● | – | – | – |
| Solid advection | – | ● axial | ● angular | ● loop | – | – | – |
| Wall energy | ● | opt | opt | opt | opt | opt | – |
| 1D liquid flow | – | – | – | – | – | ● | ● pump |
| Interface cross term | – | – | – | ● bubble ↔ emulsion | ● permeation | ● gas–liquid film | – |
| Reactions | opt | opt | opt | opt | opt | ● amine kinetics | – |
| Ports and streams | ● | ● | ● sectors | ● | ● | ● | ● |

- **Rotating bed:** a 2D mesh (axial × angle). The solid moves through the angle using the moving-bed advection block, and each sector boundary is a port. This avoids time-switching the boundary conditions.
- **Two-channel devices:** the membrane, the absorber and a heat exchanger all use the tube/shell pattern from the membrane branch. One helper covers them: `counter_current_pair(sysA, sysB, interface_cts; ncells)`.

## 5. Coupling units: ports and streams

Units with holdup (beds, tanks, nodes, boundaries) own the pressure. Flow elements (valves, pumps, compressors, vacuum pumps, pipes) own the molar flow and solve a device law. Each connection between a port and a device is a pair of `AdditiveCrossTerm`s. This is the same pattern JutulDarcy uses to couple wells to the reservoir.

Example: the current single-bed, four-stage VSA as a flowsheet.

```mermaid
flowchart LR
    Feed["Feed\nP, T, y(t)"] --> V1["Valve V1\nF, P_in, P_out"]
    V1 --> Bed["Fixed bed\nP, y, q, T, T_w per cell"]
    Bed --> V2["Valve V2\nF, P_in, P_out"]
    V2 --> Raff["Raffinate\nproduct"]
    Bed --> V3["Valve V3\nF, P_in, P_out"]
    V3 --> Pump["Vacuum pump\nF, power"]
    Pump --> CO2["CO2 product\nextract"]
```

Each stage is a schedule of valve openings on V1–V3 and set-points on the boundaries, so nothing inside the bed changes between stages. `F` can change sign; composition and enthalpy are upwinded.

- **Ports.** Each unit declares named ports with `ports(model)`. A fixed bed has `(bottom = Port(cells=[1]), top = Port(cells=[nc]))`.
- **Device laws.** A flow element has `F`, `P_in` and `P_out` as unknowns.
  - Valve: `F = Cv·x(t)·g(ΔP)`, with `g` smoothed near zero.
  - Compressor: isentropic efficiency, with power as an output variable.
- **The two cross terms in each connection.** One adds the flow (±F·y and ±F·h) to the unit's mass and energy equations. The other feeds the unit's port-cell pressure into the device's `P_in − P_port = 0` equation.
- **Stream types.** Gas, liquid and solid streams. Solid streams carry loading and temperature, which moving beds and circulating fluidised beds need.
- **Cycles become schedules.** Valve openings and feed or vacuum set-points are forces on the device models, not boundary conditions on the bed. Two-bed systems with pressure equalisation, which can't be built today, become a matter of wiring.
- **Single-unit mode stays.** The existing BCs move to `units/fixed_bed/legacy_bcs.jl` as a `PrescribedBoundary` mode, so the current examples and `test/regression/` keep working. A `PrescribedPressureSource` with the exponential ramp reproduces them exactly in flowsheet mode.

## 6. Differentiability rules

These rules are enforced in review and in tests.

1. **Anything tunable is a Jutul parameter** on a `Unit()` entity (renamed from `Column()`), or a force that can be vectorised. Nothing is read from `data_domain` inside residuals, boundary conditions or cross terms.
2. **Coefficients inside physics objects**, such as isotherm constants, are differentiated through setup with `DictOptimization`. They stay generic in `RealT` and are never typed as `Float64`.
3. **No non-smooth operations on state**, such as `max`, `abs` or `if` on pressure or flow. Use smoothed versions. Upwinding is the one accepted exception.
4. **Forces implement `vectorize_force!` and `devectorize_force`**, which makes schedules and set-points optimisable. `MultiModel` already supports `vectorize_forces`.
5. **Differentiable cyclic steady state (CSS).** Treat CSS as a fixed point `x₀ = Φ(x₀, p)`.
   - Solve it by successive substitution with Anderson acceleration.
   - Get the gradient with the implicit function theorem: an adjoint over one cycle, then a linear solve with `(I − ∂Φ/∂x₀)ᵀ`.
   - To check first: this needs the sensitivity with respect to `state0` from Jutul's adjoint.
6. **Gradient tests.** Every block and every cross term has a test comparing the AD gradient with a finite-difference gradient.

## 7. CCS-specific scope

The target is the carbon capture and storage chain: capture, then conditioning and compression. That adds:

- a component database covering CO₂, N₂, H₂O, O₂, Ar and SO₂;
- water co-adsorption for humid flue gas and direct air capture;
- a real-gas EOS for the compression and transport steps;
- process metrics in `process/`, computed from port streams and device power: purity, recovery, productivity, and specific energy in kWh per tonne of CO₂.

Specific energy needs pump and compressor models even for a pure VSA study, so equipment is worth adding early.

## 8. Migration phases

| Phase | Work | Done when |
|---|---|---|
| 0 | Merge `MoccaSystem` from the membrane branch. Rename `Column` to `Unit`. Move `data_domain` reads into parameters. Make gas density a secondary variable. | Regression tests pass |
| 1 | Generic conservation-law skeleton with source-term tuples. Move the fixed bed onto it. Add thermal and momentum options. | Regression results match to tolerance |
| 2 | Ports, `StreamCT`, and Feed, Product, Valve and VacuumPump models. Build the four-stage VSA as a flowsheet. | Flowsheet VSA matches the legacy BCs |
| 3 | Rebase the membrane onto the framework using `counter_current_pair`. Add a two-bed VSA with pressure equalisation. Port the balance functions. | Membrane tests pass; mass and energy balance closes |
| 4 | CSS solver, differentiable CSS, metrics and objectives. Update `examples/optimization.jl`. | CSS gradient matches finite differences |
| 5 | Moving, rotating and fluidised beds, solvent absorber and stripper, compressor, Peng–Robinson EOS. | One example and one gradient test per unit |

## 9. Open decisions

**Pressure–flow network or prescribed boundaries.**
Recommended: the pressure–flow network in section 5. It costs a few extra unknowns per device and needs smoothed valve laws, but it is the only option that allows units to be connected.

**How to make physics-object coefficients differentiable.**
Recommended: differentiate through setup with `DictOptimization`, which is simpler. Promoting the coefficients to Jutul parameters gives full adjoint coverage but makes the interface more complex.

**Fluidised-bed fidelity.**
A 0D well-mixed model is much less work than a 1D bubble/emulsion model. The right choice depends on which studies the fluidised bed has to support.

## Basis for this proposal

- Mocca.jl `main` at v0.1.0 (commit `70f6cd6`), and the feature branches named above.
- Jutul APIs checked in the installed source (0.4.31): `MultiModel`, `AdditiveCrossTerm`, `CTSkewSymmetry`, `add_cross_term!`, `vectorize_forces` for `MultiModel`, `solve_adjoint_sensitivities`, `WrappedGlobalObjective`, `CompositeSystem`.
- Not yet checked: support for sensitivities with respect to `state0`, which the CSS gradient in section 6 depends on.

## 10. Implementation status

Implemented on branch `worktree-mocca-architecture-doc`, one commit per phase. The original tests still pass, except for the regression references updated with issue 3.

| Phase | Status | What exists |
|---|---|---|
| 0 | Done | `MoccaSystem` hierarchy in `src/core/types.jl`; `Unit` entity (`Column` alias); `FixedBed` (`AdsorptionSystem` alias); boundary conditions read `Permeability`, `BedCrossSectionArea` and unit parameters instead of `data_domain`. |
| 1 | Done, except Ergun | One residual per conservation law with source-term tuples (`src/core/conservation.jl`); blocks in `src/blocks/`; thermal options `WithWall`, `Adiabatic`, `Isothermal(T)`. |
| 2 | Done | `FlowDevice` with `Closed`, `LinearValve`, `VolumetricFlow`; `PortStateCT` and `StreamCT`; `Flowsheet`, `connect!`, `set_boundary!`, `setup_flowsheet_model`; `Stage` and `setup_schedule`; `four_stage_vsa_flowsheet`; `examples/flowsheet_vsa.jl`. |
| 3 | Partly | `two_bed_vsa_flowsheet` with pressure equalisation. The wet flue gas processes of Krishnamurthy et al. (2014) (`src/process/wet_flue_gas.jl`): `lpp_vsa_flowsheet`, a single 13X bed with light product pressurisation, and `dual_adsorbent_vsa_flowsheet`, a silica gel bed feeding a 13X bed, with three components (CO2, N2, H2O), the paper's isotherms and Table S2 data, and examples. `parallel_stages` combines per-bed step sequences of different lengths; `setup_schedule(...; first_dt)` starts each stage with short steps. Not done: the membrane rebase and `counter_current_pair`. |
| 4 | Done, except updating `examples/optimization.jl` | Metrics from device streams (`stream_totals`, `purity`, `recovery`, `productivity`, `vacuum_pump_energy`, `specific_energy_kwh_per_tonne`); `simulate_to_cyclic_steady_state` with Anderson acceleration; `simulate_until_mass_balance`, the CO2-balance criterion of Haghpanah et al.; adjoint gradients through flowsheets, tested against finite differences. In `src/process/gradients.jl`: force vectorisation for device laws and port conditions, and `force_gradients` for the gradient with respect to each stage's device settings; `newton_cyclic_steady_state`, Newton's method on the cycle map with its Jacobian from adjoints; `cyclic_steady_state_gradient`, the gradient at cyclic steady state by the implicit function theorem (GMRES on `(I − ∂Φ/∂x)ᵀ μ = ∂J/∂x`), with respect to parameters and stage settings; adjoint objectives `StreamObjective` and `VacuumPumpObjective`; `examples/flowsheet_vsa_gradients.jl`. All are tested against finite differences (`test/gradients.jl`). |
| 5 | Not started | |

### Departures from the proposal

- **Boundaries are port conditions, not models.** An unconnected device port is a boundary with a `PortCondition` (pressure, which may follow an `ExponentialRamp`, temperature and composition). This avoids a model with no unknowns of its own.
- **A device is a one-cell model.** It holds its flow and a copy of the state at each port. Jutul forces only carry derivatives with respect to the entity they act on, so keeping every device variable in one cell is what lets the device law change per stage through forces and keep exact derivatives.
- **Cross-term rules.** Jutul adds a cross term's residual through its derivatives with respect to the target, and finds what a cross term depends on by tracing it once, at the current state. So cross terms must not hide a dependency behind a value-dependent branch (use `ifelse`), and a term that does not otherwise depend on the target touches it with zero weight (`_target_anchor`).
- **Adjoints need stored sub-steps.** Run the forward simulation with `output_substates = true` before calling `Jutul.solve_adjoint_sensitivities`; otherwise a step the solver split is linearised as one step.
- **Mocca drives the force adjoints itself.** `Jutul.solve_adjoint_forces` gave wrong gradients for flowsheets (issue 9). `force_gradients` vectorises the forces of each stage, builds the case for `Jutul.AdjointsDI.solve_adjoint_generic` itself, and passes it every sub-step the solver took as a step of its own, which also avoids issue 7.
- **Cyclic steady state gradients need a smooth cycle map.** With the abrupt start of a stage taken in one 1 s step, the first Newton solve of the stage needs about 15 iterations, and a start state changed by 1e-6 can converge to a solution 1e-3 away, so the cycle map jumps. Newton's method on the cycle map then stalls and finite differences are meaningless. Starting each stage with short steps (`first_dt = 0.05`) removes the jumps; it also halves the number of Anderson cycles for the four-stage VSA.
- **Newton's method for cyclic steady state.** It was not in the proposal. The slowest mode of the four-stage VSA cycle decays by only about 1% per cycle (spectral radius of `∂Φ/∂x` 0.99), so successive substitution needs thousands of cycles and Anderson acceleration stagnates near 1e-5. Newton converges from 1e-4 to 1e-12 in two or three iterations. Its Jacobian costs one adjoint per bed variable, which is affordable for a few hundred variables. The gradient needs a converged steady state: an error `e` in the steady state gives a gradient error of about `e/(1 − ρ)`, 100e here.

### Issues found in the existing code

Each is left at its old behaviour unless noted, so earlier results are reproduced.

1. **Heat of sorption and adsorbed-phase heat capacity only include the first component.** The loop in the old `columnflux.jl` ran over the length-1 energy equation buffer. `FixedBed(...; sorption_heat_all_components = true)` includes every component.
2. **Fixed: `PressurisationBC` had the sign of its flow reversed** relative to the other boundary conditions, and always used the feed composition, even when gas flowed out through the inlet. Its large half-cell conductance held the inlet at the set-point whichever the sign, so pressures were right, but at the start of a cycle, when the bed is above the start of the pressure ramp, it injected feed gas while the bed was venting. It now uses the same sign as the other boundary conditions and the upstream composition and temperature, as the flowsheet ports do. This had to be fixed together with issue 3: with each component conserved, removing gas of feed composition from a bed with no CO₂ drives CO₂ negative and the solver fails.
3. **Fixed: the inlet boundary conditions over-supplied components.** `AdsorptionBC` and `PressurisationBC` added `F·(y_feed − y)` on top of the inflow `F·y_feed`, which kept the total flow but not the flow of each component. Each component now enters at `F·y_feed`, and the flowsheet `legacy_inlet` option, which reproduced the old term, is removed. The regression references were updated: in the cyclic VSA, CO₂ held on the sorbent at the end falls by 3% and the average peak temperature over the last cycle by 0.6 K; in the breakthrough case, breakthrough is 0.5% later. The flowsheet and boundary-condition versions of the cyclic VSA now agree to within the time-stepping error. `examples/optimization.jl` counted pressurisation inflow with the adsorption sign, so its CO₂ recovery also changes.
4. **Adjoint gradients with respect to unit parameters fail with the legacy boundary conditions** ("cannot determine ordering of Dual tags"), on `main` as well. Forces mix unit-parameter derivatives into cell equations. Flowsheets do not have this problem.
5. **Fixed: adjoint gradients with respect to `SolidVolume` and `FluidVolume` were wrong** (about 2× and 2% off), on `main` as well. To differentiate with respect to parameters, Jutul turns the primary variables into parameters, and it stores parameters once for both time levels, so the time difference of a primary variable vanished. The loading equation, sorption uptake and heat of sorption differenced the loading directly, and the energy equation differenced the pressure directly. They now use the secondary copies `AdsorbedLoading` and `ConservedPressure`, as temperature already did through `ColumnConservedEnergy`. Forward results are unchanged. Rule: take time differences of secondary variables only.
6. **`setup_forces` drops the remainder of a stage** when its duration is not a multiple of `max_dt` (15 s stages with `max_dt = 2` become 14 s).
7. **`DictOptimization` with the default `deps = :case` gives gradients about 1% off when the solver splits the first report step.** This is a bug in Jutul's `AdjointsDI`: `evaluate_residual_and_jacobian_for_state_pair` resets `state0` to the initial state for every sub-step of report step 1, not only the first. Checking `step_info[:substep_global] == 1` instead of `step_info[:step] == 1` fixes it (column-length gradient 1.2% off → 5e-11). A Jutul-only reproduction and issue text are in `docs/upstream/`. Until it is fixed upstream, use `deps = :parameters` when the setup function only changes parameters.
8. **Fixed: the heat of sorption of every component after the first was divided by the first component's saturation capacity** in `DualSiteLangmuir`, instead of its own. It only matters with `sorption_heat_all_components = true`; for water on silica gel it gave about −800 kJ/mol instead of about −40. The first component's value is unchanged.
9. **Worked around: `Jutul.solve_adjoint_forces` gives wrong gradients for flowsheets.** These are Jutul bugs, in 0.4.29:
   - The generic `vectorize_forces!` and `devectorize_forces` store a length for every force of a model but look them up by a counter that skips absent forces, and a `MultiModel` passes each submodel the whole vector. Devices, whose port conditions are often absent, have their own methods (`src/units/equipment/device.jl`).
   - With sub-steps, the forces of a sub-step are looked up by its report step in a map indexed by sub-step, so after the first split step every stage gets the forces of an earlier one.
   - Even without sub-steps, the gradients are wrong: for the adsorption feed rate in the 5 s-stage test case it gives −3.6 against 384 from finite differences. `solve_adjoint_generic!` replaces the forces of every step with those of the case its setup function returns for the first step, which the force path builds for a single step.

   `force_gradients` avoids all three. These still need reporting upstream.

### Next steps

- Rebase the membrane branch onto `FlowChannel`, with `counter_current_pair` for the tube and shell sides.
- Step durations as design variables: durations are time steps, not forces, so they need a time scaling per stage.
- Gradients of global objectives (ratios such as purity and recovery) directly, rather than by the quotient rule from sum objectives.
- Report the force adjoint bugs (issue 9) to Jutul, with a Jutul-only reproduction.
- Ergun momentum, a Peng–Robinson EOS, and the Phase 5 units.
- Move plotting into a Makie package extension (a breaking change for `plot_*` users).

