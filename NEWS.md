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
  temperature (`pump_heat_capacity_ratio`), instead of 1.4. The presets and
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
