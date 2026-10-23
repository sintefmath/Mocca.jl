# Mocca.jl release notes

## Unreleased

### New

- **Axisymmetric beds.** A fixed bed can be meshed in the (z, r) plane of its
  column: `column_mesh(ncells, L, r_in; radial_cells, radii)` and
  `mocca_domain(mesh; axisymmetric = true, radial_diffusion_coefficient,
  radial_thermal_conductivity, ...)` give cells that are rings of the
  cylinder, with their volumes and face areas. The gas disperses and conducts
  across the radius with its own coefficients (the `RadialDispersion`
  parameter and the radial conductivity), and the wall touches the outer ring
  only (`WallCells`, `wall_cells`); inside it the wall temperature follows the
  gas. The ports are the end rings, so the bed connects to devices as a
  column does, and a bed without radial gradients gives the column's results.
  `add_bed!(fs, name, adsorbent, column; radial_cells, radii,
  radial_dispersion)` builds one, as do the presets with `radial_cells`; the
  column then gives its `radial_conductivity` and `wall_heat_transfer`, the
  coefficient between the gas next to the wall and the wall, which takes the
  place of `inner_heat_transfer`. The radial dispersion is
  `0.7 D_m + 0.1 v d_p` (`compute_radial_dispersion`, `bed_radial_dispersion`)
  at the same velocity as the axial one, unless given. Process files read and
  write these settings. `cyclic_steady_state` starts such a process from its
  steady state on one radial cell (and coarser grids along the column, as for
  any bed), spread over the rings; the bed of one radial cell takes the wall
  coefficient `1/(1/h_w + r/(4λ_r))` that gives the radial bed's mean
  temperature. `write_results` adds `bed_fields.csv`, every cell of such beds
  with its ring and radius, and the summary their radial heat transfer and
  the largest temperature difference across the radius over the last cycle.
  `plot_bed_section` draws a cut through a bed along its axis, with the wall,
  and `plot_radial_profiles` profiles across the radius. `plot_cell` and
  `plot_cell_comparison` take a port (`:top`, `:bottom`), the mean over its
  cross-section, or a cell `(i, j)` along the column and across the radius.
  The wall exchanges heat with the outer ring through the film and the half
  of the ring between its centre and the wall (`RingWallResistance`), which
  makes the cooling of a bed by conduction alone match the series solution to
  second order in the width of the rings. The manual describes the model, its
  discretisation and its verification, and the example "A column with a
  cooled wall, across its radius" compares such a bed with a 1D one.
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

- **Ports of several cells.** A port of a distributed unit is now the vector
  of cells at that end, the first or last of its `axial_layers` (one cell for
  a column), so that a unit meshed across its cross-section connects to
  devices as a column does. A stream through a port of several cells is shared
  between them by `port_weights` (their share of the cross-section); out of
  the unit each cell gives up its own gas, and the device sees the cells mixed,
  at the mixing-cup temperature, so every component and the enthalpy are still
  conserved. `port_conductance` takes the cells of a port or its name.
  Results and plots along a bed show cross-section averages (`axial_profile`,
  `axial_positions`, `cross_section_mean`, `port_state`). The face
  transmissibilities are computed from the data domain rather than the mesh,
  and conduction along the wall from its own face parameter,
  `WallTransmissibilities`. A column gives the same results as before, to
  round-off.
- **Parameters read in equations.** An equation that checks whether its
  unit has a parameter, such as the `TimeScale` of the stage-duration
  gradients, now finds it in the states Jutul evaluates equations on, whose
  fields are not the state's variables. Before, the time scale was missed
  when the equations were assembled with their derivatives; it is 1 in a
  simulation, so results do not change.
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
