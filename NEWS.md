# Mocca.jl release notes

## Unreleased

### New

- **Gas thermodynamics.** The heat capacities and the equation of state of a
  gas come from an `AbstractGasThermodynamics` model in `Mocca.Physics`:
  `IdealGasConstantCp` (the default, from `heat_capacities` per mass as
  before) or `IdealGasShomate`, with heat capacities that change with
  temperature, chosen with `Gas(names; thermodynamics = :shomate)` or by
  giving `shomate` coefficients. The common components carry the NIST Shomate
  coefficients. A new model needs one method, `component_heat_capacities`.
  A bed's system holds the model of the gas (`thermodynamics`) and of the
  adsorbed phase (`adsorbed_thermodynamics`) in place of `heat_capacity_gas`
  and `heat_capacity_adsorbed`, which `FixedBed` still takes.
- **Transport properties that follow the state.** `Gas(...; transport =
  :chapman_enskog)` lets the viscosity and the molecular diffusivity of each
  component change from their values at the feed with the temperature,
  pressure and composition of each cell, by the kinetic theory of Chapman and
  Enskog with the collision integrals of Neufeld and the mixing rules of Wilke,
  as in the thermodynamics library of the MATLAB model. The axial dispersion
  then follows the diffusivities and the gas velocity at each face, by Wakao
  and Funazkri, `ε D_L,i = 20 D_m,i + 0.5 |u| d_p`, times a `DispersionScale`
  parameter; the linear driving force takes the diffusivities at the feed
  pressure. The common components carry their Lennard-Jones parameters (all but
  CO and He). The given viscosity and diffusivities, and `FluidViscosity`,
  stay the values at the feed, so design variables and gradients act on them
  as before; a new model is an `AbstractGasTransport` with two functions. The
  default, `ConstantTransport`, gives the same results to the bit.
- **Cubic heat capacities.** `Gas(...; thermodynamics = :polynomial)`, with the
  coefficients of Reid, Prausnitz and Poling (`IdealGasPolynomial`), as in the
  MATLAB model; within 2% of the NIST values from 300 to 500 K.
- **Gas properties page.** A page of the manual's library on the gas property
  models: how to choose them, the component data and its sources, their
  accuracy against NIST, their limits and their effect on the presets.
- **Moving beds.** `Mocca.MovingBeds` adds `MovingBed`, a column whose
  sorbent moves down through it at the `SorbentMassFlow` parameter, carrying
  its loading, its heat and the gas in its pores, while gas flows through it
  either way. Its thermal
  models are `Adiabatic`, `Isothermal`, `WithWall` and the new
  `HeatExchanger`: tubes in each cell carrying a heat transfer fluid, with the
  tube wall and fluid temperatures as unknowns. `SorbentConnection`
  (`connect_solids!`) passes sorbent from one unit to the next and
  `HeatTransferFluidConnection` (`connect_heat_transfer_fluid!`) the heat
  transfer fluid, closed loops included. Two device laws come with them:
  `FixedMolarFlow` and `BackPressureController`, a blower holding the
  pressure at its inlet. All new parameters get adjoint gradients.
- **Moving beds in the flowsheet interface.** `add_moving_bed!` adds a
  section of a moving bed to a `Flowsheet`, with optional `tube_bundle` heat
  exchanger and structured packing. A `SteadyProcess` runs a flowsheet whose
  devices keep their settings, as moving beds do, to steady state
  (`steady_state`), and `performance` gives its purity, recovery, power, heat
  duties, sorbent inventory and circulation time. Process files read and
  write moving beds, `sorbent_connections`, `fluid_connections` and a steady
  `"operation"` in place of a cycle. A `Process` refuses moving beds.
- **Virial isotherm.** `Virial`, the multicomponent Virial isotherm of
  Mondino et al. (2022), and `virial` for its data by component, readable
  from adsorbent files. The library has the activated carbon of that paper,
  `activated_carbon_mondino_2022`.
- **Temperature in mass transfer, and Knudsen diffusion.**
  `compute_mass_transfer_rate` now takes the temperature, as
  `compute_mass_transfer_rate(mt, C, q, qstar, T)`; models that do not use
  it keep the four-argument method, which the new one falls back to. A bed
  whose transport properties follow the state passes the diffusivity factors
  after the temperature, `compute_mass_transfer_rate(mt, C, q, qstar, T,
  diffusivity_factors)`.
  `LinearDrivingForce` takes an optional `macropore_radius` and the
  components' `molar_masses`, which add Knudsen diffusion in series with
  molecular diffusion (`pore_diffusivity`). An `Adsorbent`, and an adsorbent
  file, may give a `macropore_radius`.
- **Limited second-order advection in moving beds.** What the sorbent
  carries (loading, heat, pore gas) and the heat transfer fluid in the tubes
  are reconstructed at the cell faces with a slope limiter, `VanAlbada()` by
  default, set by the `limiter` keyword of `MovingBed` and `add_moving_bed!`
  and in process files. `Upwind()` is the former first-order scheme, and
  `VanLeer()` and `Minmod()` are also available; a new limiter is a subtype
  of `AbstractLimiter` with a method of `limited_slope`. First-order
  upwinding added a numerical axial conduction of about 10 W/(m·K) to the
  MBTSA cooling section at 400 cells per metre, 50 times the physical one,
  and the plant's capture rate drifted with the grid; with the limiter it
  is grid-converged from 100 cells per metre.
- **Metrics of moving bed processes.** `stream_rates`, `heat_duty`,
  `sorbent_inventory`, `sorbent_circulation_time` and `sorbent_path`.
- **Moving bed TSA of Mondino et al. (2022).** The "Moving bed TSA" example
  builds the five-section MBTSA process for the flue gas of a waste-to-energy
  plant, and `plot_moving_bed_profiles` draws the profiles through the
  sections. With the purge at the paper's 15% of the flue gas and a
  water-side heat transfer coefficient of 10000 W/(m²·K), which the paper does
  not give and which matches its heat loads, it gives 99.7% purity, 93.1%
  capture and 5.53 MJ/kg, against the paper's 97.2%, 90.8% and 5.7 MJ/kg, with
  heat loads within 1.5%. The capture rate hangs on a threshold in the purge:
  below about 7.9 kg/s no purge gas is left to strip CO₂ from the hot sorbent
  in precooling, and it falls to 87%.

- **Gradient checks.** `check_gradients(p)` compares the gradients of a
  process with central differences of its steady state at two step sizes,
  and marks each pair of variable and measure as passing, inconclusive (the
  two differences disagree, as at a kink or a jump) or failing (they agree
  with each other but not with the gradient). `test/gradient_checks.jl` runs
  it on every process of the library and every process the examples build,
  with `MOCCA_GRADIENT_CHECKS=true`. The extending guide has the rules for
  differentiable models.
- **Optimisers.** `optimise!` meets its targets on measures with an
  augmented Lagrangian (`constraints = :lagrangian`, the default): after each
  BFGS run the multiplier of a target grows by the shortfall left, and the
  penalty tenfold where the shortfall did not shrink, until every target is
  met to `constraint_tol`. Its first run is the former penalty method,
  `constraints = :penalty`, which stops there and gives the former results.
  Each design is evaluated once (the restarts reuse it). `method = :pattern`
  searches without gradients, by a compass search in the box, for
  objectives with kinks or jumps. A design where the process does not settle
  or fails to simulate is flagged in the history (`settled`), left out of the
  best, and the optimiser steps back from it. The result has `feasible`,
  `multipliers` and, per evaluation, `settled` and `violation`.
  `pareto_front!` traces a Pareto front by ε-constraints, each optimum
  started from the last, and `explore` evaluates designs spread over the box
  (a Halton sequence), each from the initial state, to find the regimes and
  a start before following gradients. A `SteadyProcess` starts each
  `steady_state` from the last one it found (`warm_start = true`).
  In the "Optimization" example the targets are now met exactly: 65.4%
  recovery at 90.0% purity (65.5% at 89.9% before), and 115 kWh/t at 90.0%
  purity and 60.0% recovery (116 kWh/t at 89.9% and 59.9% before), in 28 and
  56 evaluations (20 and 24 before).
  Neither is the end of the flat ridge the recovery rises along: one BFGS
  run of 60 iterations reaches 68.0%, with the adsorption time at its bound.
  The library has the MBTSA, `mondino_2022_mbtsa`, at 50 cells per metre,
  and a new example, "Optimising the moving bed TSA", finds its least
  equivalent work for 95% purity and 90% capture over the purge, the hot
  water temperature, the sorbent flow and the height of the desorption
  section, and its trade-off with productivity.
- **Gradients for any place, and geometry as a design variable.** The
  change of the model's inputs with a design variable comes from setting it
  a little up and down, without simulating: the stage settings and
  durations, now also every unit parameter, and the quantities such as the
  adsorbent volume, by central differences with a power-of-two step, exact
  for the built-in places, whose effect is linear (the gradients of step and
  feed variables move by about 4e-11 as the rounding of the former step
  goes). The adjoints then take only the parameters that change. A new kind
  of place needs only a way to read and write its value. New places:
  `GeometryPlace(unit, field)`, a field of a bed's or moving bed section's
  `Column`, such as its length, radius or porosity, or the packing fraction,
  which rebuilds the unit on as many cells (a new radius keeps the wall
  thickness); `ScaledPlace(place, factor)`, a place set to a multiple of the
  variable; and `SettingPlace(device, field; port)` for the pressure or
  temperature at a boundary port, such as the extraction pressure. `Place`,
  their abstract type, is exported and documented with what a new kind
  implements. Gradients
  match central differences to within 1e-5 for the column length and radius
  of the four-stage VSA and the height of the MBTSA's desorption section and
  its extraction pressure. Process files read and write the new places, and
  the variables of steady processes.
- **Gradients of steady processes.** `performance_gradient` and `optimise!`
  work on a `SteadyProcess`, such as the MBTSA. A steady state is the fixed
  point of a cycle of one long backward-Euler step, which returns any state
  that satisfies the steady equations, so `cyclic_steady_state_gradient`
  gives its gradients unchanged: the adjoint of the steady equations at the
  steady state, linearised there without solving the step (its new `steps`
  keyword). Variables set unit parameters (`ParameterPlace`), such as the
  sorbent flow of every section or the water temperatures, or device
  settings, with the new `SettingPlace(device, field)` for a field of a
  device's flow law, such as the purge rate. On the MBTSA at 10 cells per
  metre the gradients of six measures with respect to the sorbent flow, the
  hot and cold water temperatures and the purge match central differences
  of the steady state to within 4e-6, and at the example's 200 cells per
  metre within 2e-5. The MBTSA example gains a section on the sensitivities
  of its capture rate and equivalent work to the sorbent flow, the hot water
  temperature and the purge. `steady_state` takes `progress`.
- **Measures as functions of integrals.** Every performance measure is now
  a function of integrals over the cycle (rates at steady state), such as
  the moles of the key component in the extract or the pump work, and of
  quantities that do not depend on the state, such as the cycle time:
  `performance_integrals`, `performance_quantities`, `performance_measures`
  and `performance_sums` list them. `performance_gradient` takes the chain
  rule from the integrals to each measure by automatic differentiation
  (`measure_partials`), in place of the quotient rules written out by hand,
  and gives the same gradients to rounding. `add_measure!` adds a measure of
  the user's own, from these or from an integral added with `add_integral!`
  (any Jutul sum objective), which `performance` reports and
  `performance_gradient` and `optimise!` differentiate with no derivative
  code. `performance_gradient` takes `measures`, the measures to
  differentiate. When
  the measures depend on more integrals than there are measures, each
  measure takes one adjoint of its linearisation instead of one per integral.
- **Fan and blower work on links.** A `:link` device that moves gas up in
  pressure, such as the recycle and purge blowers of the MBTSA, counts as a
  fan or blower: `BlowerObjective` and `blower_energy`, with the pumps'
  efficiency and heat capacity ratio. Its work is part of the `energy` (per
  cycle) and `power` (steady), and reported as `blower_energy` and
  `blower_power`. A link that lets gas down, such as the valve of an
  equalisation, takes none, so VSA results do not change.
- **Heat duty and its equivalent work as integrals.** `HeatDutyObjective`
  gives the heat `ṁ·c_p·(T_in − T_out)` of a unit's heat transfer fluid for
  adjoints, and `EquivalentWorkObjective` (`equivalent_work`) the work it is
  worth against the new `ambient_temperature` option (298.15 K) of `Process`
  and `SteadyProcess`: each bit of heat weighted by the Carnot factor at the
  temperature the fluid gives it up at, `ṁ·c_p·[(T_in − T_out) − T₀·ln(T_in/T_out)]`.
  The performance of a process with heat exchangers adds
  `equivalent_work` and `specific_equivalent_work`, the electrical and
  equivalent work together per tonne of key component. `heating` and
  `cooling` now count the units whose fluid enters from outside above or at
  and below the ambient temperature, which for the MBTSA is the same as
  before. Objectives now see each unit's parameters in the state in every
  adjoint (`objective_state`), as Jutul's own adjoints give them.
- **`write_results`.** Writes a process result to CSV files for a
  spreadsheet, in bara, °C, mol%, kg/h, Nm³/h and mol/kg: a summary with the
  performance and the mass balance, a stream table per cycle and per stage,
  the streams and bed ends through time, and profiles along each bed at the
  end of each stage.
- **Progress of the cyclic steady state.** `cyclic_steady_state` prints each
  cycle as it runs: the cycle count, the change over the cycle against the
  tolerance, a bar of how far it has come, and an estimate of the cycles and
  time left. Each phase has its own line and bar: the plain cycles and then
  Newton's method, with each Jacobian or GMRES solve as it goes and each grid
  of a grid sequence, or the mass balance. In a terminal the lines are drawn
  together and rewritten in place. It is on by
  default in an interactive session and set with `progress`, which the
  solvers `simulate_to_cyclic_steady_state`, `newton_cyclic_steady_state` and
  `simulate_until_mass_balance` also take. Newton's method no longer prints
  its iterations at `info_level ≥ 0`.

### Changed

- **Diffusivities at the feed.** `mass_transfer(adsorbent, gas, feed)` and
  `axial_diffusivity(gas, feed)` take the feed, at which a transport model
  gives the diffusivities. Without it they still work for a gas with given
  or constant diffusivities, and raise an error for one whose diffusivities
  come from a transport model. A bed whose dispersion follows the state, by
  Wakao and Funazkri, has a `DispersionScale` parameter in place of
  `AxialDispersion`, and `bed_dispersion` raises an error for it.
- **Molar heat capacities in the energy balance.** The heat capacity of the
  gas is the mole average of the molar heat capacities, `Σ yᵢ Mᵢ cₚᵢ`. It
  was `(Σ yᵢ cₚᵢ)·(Σ yᵢ Mᵢ)`, exact only for one component or for equal
  heat capacities per mass: for CO2/CH4 at 40/60 it gave 45.6 J/(mol·K)
  instead of 36.3. The adsorbed phase is weighted by its loadings rather than
  by the gas composition: a heat capacity that followed the gas over the
  cycle added heat on every cycle, and the wet flue gas presets of
  Krishnamurthy et al. (2014) ran 12 K above the feed temperature on average
  at the cyclic steady state; they now run at it (LPP VSA recovery 78.3% →
  78.8%, dual-adsorbent VSA 69.8% → 65.5%, specific energy 188 → 196 kWh/t). The sorption term and the streams into a bed give
  back exactly the `C_pg·T` the gas terms count for the moles they move, so
  the energy balance stays right when heat capacities differ between
  components or change with temperature; a stream entering a bed brings the
  integral of the heat capacity from the bed's temperature to its own. With
  the constants of Haghpanah et al. (2013), 30.7 J/(mol·K) for both
  components, the gas heat capacity is now that constant, where it was up to
  5% higher in a mixture.
- **Heat capacity ratio of the pumps from the feed.** A `Process` without
  `heat_capacity_ratio` takes `γ = c_p/(c_p − R)` of the feed at its
  temperature (`pump_heat_capacity_ratio`), instead of 1.4, for its vacuum
  pumps, compressors and blowers, cyclic or steady. The presets and
  library processes set 1.4, as their papers do.
- **`TimeScale` in every residual.** Units checked for the `TimeScale`
  parameter with `hasproperty`, which is false on the local
  automatic-differentiation states Jutul assembles the equations with, so it
  was left out there; it is 1 except when differentiating with respect to a
  stage duration, so simulations are unchanged.
- **Renamed variables.** `SpecificHeatCapacityFluid` and
  `SpecificHeatCapacityAdsorbent` are now `GasMolarHeatCapacity` and
  `AdsorbedMolarHeatCapacity`; they are still `:C_pg` and `:C_pa` in the
  state.
- **`Process` is the abstract type of every process.** The cyclic process
  is now `CyclicProcess`, and its result `CyclicResult`; `SteadyProcess` is a
  `Process` too, and `SteadyResult` a `ProcessResult`. `Process(fs, cycle)`
  still builds a cyclic process, every `p::Process` annotation accepts both
  kinds, and process files do not name the type, so existing code and files
  work unchanged. Both kinds share one interface: `steady_state` finds the
  state a process settles to (the cyclic steady state of a cyclic process),
  `stages` its operation, `performance(p)` measures it there, and design
  variables, measures and integrals are added the same way; a steady process
  takes variables at `ParameterPlace`s, with gradients still to come. A
  cyclic process shows as "Cyclic process". The extending guide lists what a
  new kind of process implements.

- **Gas balances converge against the gas flowing through each cell.** The
  mass and energy errors of a cell are measured against the gas it holds plus
  the gas that flows through it over the step, and its heat capacity. Before,
  only the gas held counted, and on fine grids with a fast flow the round-off
  in the flux between two nearly equal pressures kept the error above the
  tolerance on long steps: the solver cut every 60 s step five times, and the
  moving bed TSA at 1280 cells per metre took 200 s instead of a few seconds.
- **Cyclic steady state from a few plain cycles.** `cyclic_steady_state`
  (with the default `method = :newton`) now starts Newton's method after
  `handover = 3` plain cycles, instead of after Anderson acceleration to a
  cycle change of `1e-3`, which was most of the time of a cold solve.
  `newton_cyclic_steady_state` controls its steps: a step from which the
  cycle fails to simulate, or that makes the change more than `growth` (100)
  times larger, is cut back, up to `max_backtracks` (6) times. If Newton's
  method fails or does not converge, Anderson acceleration goes on to `1e-1`
  and then `1e-3`, and Newton's method tries again from each. `handover = 1e-3`
  gives the earlier behaviour. Cold solves to `1e-9` on 10 cells: the
  four-stage VSA in 4 s instead of 26 s, the wet flue gas LPP VSA in 30 s
  instead of 437 s, and the wet flue gas dual-adsorbent VSA in 37 s, where it
  failed before.
- **Cycle Jacobian from forward sensitivities.** Newton's method builds the
  Jacobian of the cycle in one forward sweep over the cycle, every column at
  once, instead of one adjoint per variable of the beds. A Jacobian takes
  about as long as a cycle: Newton's method on the four-stage VSA takes 1.3 s
  instead of 26 s on 10 cells, and 3.3 s instead of 120 s on 30 cells.
- **GMRES on fine grids (`krylov`).** `newton_cyclic_steady_state` and
  `cyclic_steady_state` can solve each Newton step by GMRES, one forward
  sensitivity per product, without building the Jacobian. The default,
  `krylov = :auto`, does so on grids whose beds have more than 800 variables
  (a bed of about 130 cells with two components) and builds the full
  Jacobian below; `krylov = true` or `false` fixes the choice, and
  `krylov = n` moves the threshold to `n` variables.
- **Grid sequencing.** On a bed of more than `coarse_cells = 40` cells,
  `cyclic_steady_state` first finds the steady state on beds of half as many
  cells, down to 40 or fewer, and starts each finer grid from the coarser
  one's steady state, interpolated. The wet flue gas LPP VSA on 160 cells now
  converges cold in about 3 minutes, where Newton's method stalled for over
  50 minutes before. `coarse_cells = nothing` starts every grid from its
  initial state.
- **Kept time steps near the steady state.** Once the cycle change is below
  `keep_steps_below` (`1e-3`), Newton's method keeps the time steps of its
  cycle, including any the solver cut, for every later cycle, so that the
  cycle map it works with is smooth.
- **Cyclic steady state judged against what the steps resolve.** The
  change over a cycle (`cycle_change`) is now a mixed absolute and relative
  test, as in an ODE solver: an entry has settled when its change is at most
  `tol` times its scale plus what steps at the nonlinear tolerance resolve
  (that tolerance for a mole fraction, times the cell's total loading for a
  loading, times the pressure for a pressure, and in kelvin for a
  temperature). The change reported is the largest excess over that
  resolution relative to the scale, so it is below `tol` exactly when the
  test passes. Every component's mole fraction is checked, and Newton's
  method weighs its residual the same way. Before, a trace as the last gas
  component, which only moves by minus the change of the others, was
  measured against a floor of 10⁻⁶ and held Newton's method at a change of
  exactly 10⁻⁶. On 30 cells, with the gas in the order CO2, H2O, N2, the wet
  flue gas LPP VSA now converges in 21 cycles (28 s) and the dual-adsorbent
  VSA in 82 (4 min), where both stopped unconverged after 98; dry flue gas
  with a trace of water last converges in 19 instead of not at all. In the
  usual order the library processes take about as many cycles as before.
  `cycle_change` without `nonlinear_tolerance` is the plain relative change,
  as before.
- **One nonlinear tolerance for a cyclic steady state.** With
  `method = :newton`, the plain cycles and any Anderson acceleration before
  Newton's method are simulated at Newton's tolerance, `1e-8`, that of the
  adjoints, instead of `1e-6`, so that every phase works with the same cycle
  map. `cyclic_steady_state` takes `nonlinear_tolerance` to change it.
- **Inventory closure.** A `CyclicResult` has `closure`: for each component,
  the change of what the beds hold over the last cycle as a fraction of the
  moles of it fed. It is printed at the end of the progress of
  `cyclic_steady_state`, and does not decide convergence.
- The four-stage, two-bed and wet flue gas presets give their feed by
  component name, so that their `gas` may list the components in another
  order.
- `cyclic_steady_state` no longer fails with a `FieldError` when Newton's
  method fails from a hand-over that is not the last; Anderson acceleration
  goes on to the next, as documented.
- **Submodules.** The source is split into one submodule per layer
  (`Mocca.Framework`, `Mocca.Physics`, `Mocca.Balances`, `Mocca.Devices`,
  `Mocca.FixedBeds`, `Mocca.Networks`, `Mocca.Operation`, `Mocca.Interface`),
  plus `Mocca.Legacy` for 0.1 inputs and `Mocca.Plotting`. `using Mocca`
  exports the same names as before, and `Mocca.X` still works for every name
  with a docstring. Other internal names are now reached through their module,
  such as `Mocca.Operation._simulate_steps`, and names from CairoMakie, such as
  `Mocca.Figure`, are no longer reachable through `Mocca`.
- Loading Mocca no longer creates a `results` directory.
- **Examples.** Each example now covers its own part of the package:
  `getting_started` (now with the CSV export), `process_from_parts` (which
  replaces `custom_setup_cyclic_vsa` and `cyclic_vsa_haghpanah_2013_co2_n2`),
  `breakthrough` (formerly `dcb_haghpanah_2013_co2_n2`, now comparing the
  thermal models), a new `two_bed_vsa` with pressure equalisation, and
  `unit_network` (formerly `flowsheet_vsa`).

## v0.2.0

Mocca 0.2 is a rewrite of the simulator. A process is now a flowsheet of beds
joined by flow devices, run through a cycle of steps. Mocca finds its cyclic
steady state and gives adjoint gradients of its performance. This release is
breaking: see **Migrating from 0.1** in the docs.

### New

- **Flowsheets.** Beds, flow devices (closed valves, linear valves, flow
  controllers) and links between beds, coupled so that each component and the
  stream enthalpy are conserved. Multi-bed processes such as a two-bed VSA with
  pressure equalisation are supported.
- **Cycles of steps.** `Pressurise`, `Adsorb`, `Blowdown`, `Evacuate`, `Purge`,
  `Equalise` and `Idle`, run in sequence or concurrently across beds. Steps find
  their devices by role, and pressure ramps carry over from step to step. You
  can add your own kinds of step.
- **`Process`.** One object holds the flowsheet, the cycle and the solver
  options, with `simulate`, `cyclic_steady_state` and `performance`.
  Performance covers purity, recovery, specific energy and productivity.
- **Faster cyclic steady state.** Anderson acceleration, then Newton's method
  with an adjoint Jacobian.
- **Gradients and optimisation.** Adjoint gradients at cyclic steady state with
  respect to device settings, stage durations and model parameters. Expose them
  as named design variables with `add_variable!`, then use
  `performance_gradient` and `optimise!`.
- **Process description.** `Gas`, `Feed`, `Adsorbent`, `Column` and `Wall`, all
  in SI units, with conversions in `Mocca.Units` (e.g. `0.1bar`).
- **Files and library.** JSON files for processes and adsorbents, and a library
  that includes:
  - Haghpanah et al. (2013): the four-stage VSA
  - Krishnamurthy et al. (2014): the wet flue gas LPP and dual-adsorbent VSA
    processes
- **Thermal models:** `WithWall` (default), `Adiabatic` and `Isothermal`.
- **Input checking.** `validate(process)` lists likely input mistakes without
  stopping.
- **Plots** of flowsheet streams, bed profiles and comparisons.
- **Extension points** for isotherms, mass transfer models, source terms,
  units, device laws, steps and connections.

### Changed

- **Productivity** is in mol of key component per m³ of adsorbent per second.
- **Loadings** are in mol per m³ of adsorbent throughout. `equilibrium_loading`
  and `isotherm_table` report mol/kg, for comparison with published isotherms.
- **Convergence checks are scale-free:**
  - mass errors are measured relative to each cell's holdup;
  - energy errors are measured as a temperature change in K;
  - cycle-to-cycle changes are measured per component.

  The default nonlinear tolerance is 1e-6.
- **Dependencies:** Jutul 0.4.29–0.5.2 is required, and Tullio is no longer a
  dependency.

### Fixed

These fixes change results compared with 0.1:

- **Pressurisation flow sign.** The pressurisation boundary condition had its
  flow sign reversed and always used the feed composition. When the bed started
  a cycle above the boundary pressure, feed gas was pushed in instead of gas
  venting out. The flow now follows the pressure difference and carries the
  upstream composition and temperature.
- **Feed inlet.** The adsorption and pressurisation inlets added a correction
  to the feed inflow that kept the total flow right but not each component's.
  Each component now enters at its feed mole fraction.
- **Heat of sorption.** It and the adsorbed-phase heat capacity now include
  every component, not just the first. In the Haghpanah VSA, breakthrough comes
  2.2% later and the bed runs 1–2 K cooler.
- **Heat of adsorption.** Each component's value is now averaged over its own
  site capacities. Before, the first component's capacities were used for every
  component.
- **Adjoint gradients** with respect to solid and fluid volume are now correct.

### Removed

- **The 0.1 interface:**
  - the `adsorptionConstants` workflow;
  - the per-stage boundary condition types (`AdsorptionBC`,
    `PressurisationBC`, `BlowdownBC`, `EvacuationBC`);
  - `setup_boundary_conditions` and `setup_mocca_case`.
- **Migrating inputs:** `convert_constants` turns a 0.1 constants struct or
  input file into the new objects.
