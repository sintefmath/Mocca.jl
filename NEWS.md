# Mocca.jl release notes

## Unreleased

### Changed

- **Submodules.** The source is split into one submodule per layer
  (`Mocca.Framework`, `Mocca.Physics`, `Mocca.Balances`, `Mocca.Devices`,
  `Mocca.FixedBeds`, `Mocca.Networks`, `Mocca.Operation`, `Mocca.Interface`),
  plus `Mocca.Legacy` for 0.1 inputs and `Mocca.Plotting`. `using Mocca`
  exports the same names as before, and `Mocca.X` still works for every name
  with a docstring. Other internal names are now reached through their module,
  such as `Mocca.Operation._simulate_steps`, and names from CairoMakie, such as
  `Mocca.Figure`, are no longer reachable through `Mocca`.
- Loading Mocca no longer creates a `results` directory.

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
