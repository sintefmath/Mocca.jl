# A simpler user interface for Mocca

*Design note — draft for discussion. 2 October 2026. Based on branch `worktree-mocca-architecture-doc` (commit `f83b0ec`); follows [the architecture proposal](mocca-architecture-proposal.md).*

*Status, 2 October 2026: implemented in `src/interface/`, with the decisions of section 5 taken as JSON, plain SI with helpers, deprecation of the old interfaces, and an open set of steps. Section 6 lists what was built and where it departs from this note.*

The architecture proposal made Mocca a flowsheet simulator: units, devices, ports and schedules. The layers underneath work, but setting up or changing a case still takes a lot of low-level code. This note plans a layer on top that lets users add and edit components, isotherms, adsorbents, columns, flowsheets and cycles in process terms. The solver, gradients and tests stay as they are.

## 1. What makes it hard today

- **One large settings struct holds everything.** `HaghpanahConstants` (`src/io/constants.jl`) mixes isotherm, adsorbent, column, wall, gas and operating settings. It is fixed at two components, with fields such as `molecularMassOfCO2`, and the wet flue gas setup (`src/process/wet_flue_gas.jl`) uses it for three components. A second struct, `adsorptionConstants`, and `processInfo` overlap with it.
- **Isotherm parameters are hard to edit.** `DualSiteLangmuir` takes parallel arrays per adsorption site, ordered by component (`b0`, `qsb`, `ΔUb`, ...). Changing one component means editing six arrays in step, and nothing checks that they line up.
- **Physical properties are spread around.** An adsorbent's properties sit partly in `FixedBed`, partly in the constants and partly in the domain. `mocca_domain` takes 16 keywords. Users compute permeability and dispersion themselves (`examples/custom_setup_cyclic_vsa.jl`), and dispersion quietly depends on the feed velocity through `V0_inter`.
- **Flowsheets take a lot of code**, about 30 lines for one bed. Users have to:
  - repeat the component names for every device (`setup_device_model(names)`);
  - spell out pressure, temperature and composition at every boundary (`PortCondition`);
  - work out valve conductances in m³/(s·Pa) (`port_conductance`) and flow rates (ε·v·A).
- **Stages are verbose and low-level.** Each stage is nested settings per device, with pressure ramps written out each time. The same feed composition and temperature appear many times.
- **Ready-made processes are hard to modify.** `four_stage_vsa_flowsheet`, `two_bed_vsa_flowsheet`, `lpp_vsa_flowsheet` and `dual_adsorbent_vsa_flowsheet` are 60–100-line functions. To change one step, a user copies the whole function.
- **Running a case means juggling five objects:** model, parameters, initial state, forces and time steps. Metrics and objectives need device names such as `:V_vacuum`, and performance is assembled by hand in each example.
- **Design variables are mapped by hand.** In `examples/haghpanah_2013_moo/problem.jl`, `design_gradient` adds up internal slots by hand: P_I is the blowdown end plus the evacuation start, P_L the evacuation end plus the pressurisation start. Every new study would repeat this.
- **There are two front ends.** The old boundary conditions (`AdsorptionBC` and the others, with `setup_forces`) and the flowsheets both exist. JSON input (`parse_input`, `setup_mocca_case`) only supports the old one.

## 2. Principles

1. **Separate the physics from the plant from the operation.** Adsorbent data should not know about valves, and a cycle should not know about cell counts.
2. **Describe things by component and by name**, not by index or parallel array.
3. **Use defaults derived from physics**, so most settings can be left out. For example, a valve's size defaults to one that reaches the set pressure in about the ramp time.
4. **Speak the user's language:** velocities, bar, kJ/mol and step names, converted internally.
5. **Build a new layer on top, not a rewrite.** Every new object builds today's `Flowsheet`, `Stage` and `FixedBed`, so the solver, the gradients and the tests are untouched.

## 3. Layers

### 3.1 Components and gas

- A small registry of named components (CO₂, N₂, H₂O, ...) with molar mass and gas heat capacity.
- A gas mixture is an ordered list of components. It sets the component order everywhere else.
- Feeds and other compositions are given as named fractions, such as "CO₂ 0.15, N₂ 0.85", and checked to sum to 1.

### 3.2 Isotherms

- Entered per component, with one entry per site (saturation, affinity, energy) and explicit units.
- The mixture isotherm is assembled from the per-component entries in gas order. A component without an entry adsorbs nothing, which removes the zero-padding needed today.
- One common interface (today's `AbstractIsotherm`), so Langmuir, Sips, Toth and fitted tables can be added later.
- Helpers to plot and tabulate equilibrium loading against pressure and temperature. They catch unit mistakes early, such as mol/kg against mol/m³.

### 3.3 Adsorbents

One named object holding:

- the isotherm and the mass transfer model;
- density, heat capacity, particle size, particle porosity and tortuosity;
- the heat-of-sorption choice (all components or the first only);
- an optional literature reference.

The library ships 13X and silica gel (Haghpanah et al. 2013; Krishnamurthy et al. 2014) as data files. Users copy and edit a data file rather than Julia code.

### 3.4 Columns and units

- Column geometry (length, radii, bed porosity), wall and heat transfer, and the thermal model, as one object.
- A bed is then "adsorbent + column + number of cells". Permeability, dispersion and the mesh are computed inside, with an optional override for each.
- Later units (moving bed, membrane, ...) follow the same pattern: material data plus geometry plus discretisation.

### 3.5 Flowsheets

- **Devices with roles:** feed, product, vacuum pump, valve, link and equalisation. The role gives:
  - the default boundary condition;
  - which streams count as feed, product and extract for metrics;
  - a sensible default valve size.
- **Shorter construction.** Connections such as "bed bottom to the feed", without device port names where the direction is obvious. The gas is taken from the flowsheet, so devices don't repeat it.
- **Inspection.** A printed summary and a diagram of units, devices and boundaries.

### 3.6 Steps and schedules

- **Step templates in process terms:** pressurise, feed, blowdown, evacuate, purge, light-product pressurisation, equalise and idle. Each takes a target pressure, a duration and which end of the bed. It fills in the ramps, temperatures and compositions from the flowsheet, and closes every other device.
- **A cycle is a list of steps.** Multi-bed cycles use per-bed sequences combined by the existing `parallel_stages`.
- **Time stepping becomes one options object:** maximum step, first step and whether to keep the step pattern fixed. "Keep step pattern" replaces today's `reference_durations`, which only optimisers over durations need.

### 3.7 A Process object

One object holding the flowsheet, the cycle, the initial state and solver options.

- **Initial state** defaults to the bed saturated with the light product at the low pressure.
- **Simulation:** simulate for a number of cycles, find the cyclic steady state (Newton by default, the mass-balance criterion as an option), and report performance.
- **Performance** returns purity, recovery, specific energy and productivity, using the device roles. It needs no device names.

### 3.8 Design variables and optimisation

- **Named process variables** declared on the Process: P_L, P_I, feed velocity, step durations and adsorbent properties.
- **Each variable records every place it acts.** P_L, for example, is both the end of evacuation and the start of pressurisation. The chain rule over those places is then automatic, which replaces the hand mapping in the paper study.
- **Ready-made objectives and constraints** follow from this: maximise recovery subject to purity ≥ x, minimise energy, and so on.
- **One optimise call** takes bounds and wraps `Jutul.unit_box_bfgs` with the steady-state warm start. The paper study's `problem.jl` would shrink to the problem statement.

### 3.9 Input files

- One human-editable format mirroring sections 3.1–3.7: components, adsorbents, columns, flowsheet, cycle and options.
- Read and written, so a case set up in Julia can be saved.
- Library entries are files in the same format.
- It replaces the two current JSON variants and the old boundary-condition path.

### 3.10 Errors and validation

Check at construction time, not inside the solver:

- unknown component names;
- compositions that don't sum to 1;
- a step that sets a connected port;
- a valve on the wrong bed end;
- a pressure outside the cycle's range.

Messages name the unit, step and field.

## 4. Migration order

| Step | Work | Outcome |
|---|---|---|
| 1 | Components, gas, per-component isotherms, adsorbent and column objects. `HaghpanahConstants` becomes a converter into them. | Physics editable without the large struct |
| 2 | Device roles and boundary defaults, simpler connections, default valve sizing | A flowsheet in about 10 lines |
| 3 | Step templates, cycle object, time-stepping options | Cycles written in process terms |
| 4 | The Process object: simulate, steady state, performance | One object instead of five |
| 5 | Named design variables, objectives and optimise | Paper studies become short scripts |
| 6 | File format and library files; retire the old JSON and boundary conditions | One front end |
| 7 | Rewrite the ready-made processes and examples on the new layer | They become short and serve as templates |

Each step keeps the current lower-level functions working, so examples can move over one at a time. The check for each step: the four-stage VSA, two-bed VSA and wet flue gas cases give results identical to today's.

## 5. Open decisions

- **Input format.** TOML is easier to hand-edit; JSON matches the existing schema files in `models/json/`.
- **Units.** Plain SI with helpers (bar, kJ/mol), or a units package such as Unitful. Unitful catches more mistakes, but adds a dependency and may get in the way of automatic differentiation.
- **How far to deprecate.** Whether `HaghpanahConstants`, `adsorptionConstants` and the old boundary conditions (`setup_forces` and the `*BC` structs) stay as compatibility shims or are removed at the next breaking release.
- **Scope of step templates.** A fixed set of VSA/PSA steps, or an open set that users extend. The latter would also cover TSA and moving beds later.

Decided: JSON; plain SI with helpers (`Mocca.Units`); the old interfaces are deprecated shims, kept until the next breaking release; the set of steps is open.

## 6. Implementation

In `src/interface/`, one file per layer, built on the existing `Flowsheet`, `Stage`, `FixedBed` and gradient code, which are unchanged apart from new fields on `Flowsheet`.

| Step | Files | What exists |
|---|---|---|
| 1 | `units.jl`, `gas.jl`, `materials.jl` | `Mocca.Units` (`bar`, `kJ_per_mol`, `celsius`, ..., and unit strings in files); a component registry with CAS numbers (`Component`, `component`, `register_component!`); `Gas`, `composition`, `Feed`; isotherm data by component (`dual_site_langmuir`, `langmuir`, per m³ or per kg), `Adsorbent`, `Column`, `Wall`; `equilibrium_loading`, `isotherm_table`; `convert_constants` from the old structs. |
| 2 | `builder.jl` | `Flowsheet(gas, feed)`, `add_bed!` (permeability, dispersion and gas properties computed), `add_device!` with roles (feed, light inlet, product, extract, waste, link), role default boundaries and valve sizes, `print_flowsheet_diagram`. |
| 3 | `steps.jl` | `AbstractStep` and the steps `Pressurise`, `Adsorb`, `Blowdown`, `Evacuate`, `Purge`, `Equalise`, `Idle`; `Concurrent`; `Cycle` with parallel sequences; `compile_cycle`; `TimeStepping`; `register_step!` and helpers for step authors. |
| 4 | `process.jl` | `Process`, `simulate`, `cyclic_steady_state` (Newton with warm start, Anderson, or the mass balance criterion), `performance`, `initial_state`. |
| 5 | `variables.jl` | `add_variable!` with step, feed and parameter places, `p[:name]`, `performance_gradient`, `optimise!`. |
| 6 | `files.jl`, `models/library/` | `read_process`, `write_process`, `load_adsorbent`, `save_adsorbent`; library adsorbents (Haghpanah et al. 2013; Krishnamurthy et al. 2014, 13X and silica gel) and processes. |
| 7 | `presets.jl`, `examples/` | `four_stage_vsa_process`, `two_bed_vsa_process`, `lpp_vsa_process`, `dual_adsorbent_vsa_process`; the examples rewritten on the new layer, except `history_matching.jl`. |

`validate` (in `validation.jl`) and `MoccaInputError` (in `errors.jl`) give the checks of section 3.10.

Checks (`test/interface.jl`): the four ready-made processes give the same parameters, forces, initial states and simulated states, to the last bit, as the deprecated builders; every process round-trips through a file exactly; the gradients of purity, recovery, specific energy and productivity with respect to step pressures, the feed velocity and pressure, stage durations and a bed parameter match finite differences to 1e-6 or better.

### Departures from the note

- **`Column` changed meaning.** It was an alias of the `Unit` entity; it is now the column geometry. `Column()` with no arguments still returns `Unit()`, with a deprecation warning.
- **Default valve size.** A valve defaults to the conductance of the half cell at the bed end it is attached to (in series for a link), which reproduces the old boundary conditions, rather than a size that reaches the set pressure in the ramp time.
- **Pressure carries over between steps.** A step's ramp starts where the step before left its bed, around the cycle, unless `from` is given; `Equalise` leaves both beds at the mean of their start pressures. Extract and waste lines default to the lowest pressure of their bed.
- **Dispersion follows the adsorption velocity once.** With `dispersion_velocity = :adsorption` (the default) a bed's axial dispersion is set from its adsorption step's velocity when the `Process` is built, and stays fixed when design variables change, as the gas viscosity and density do. It is the bed's `AxialDispersion` parameter, so setting it does not rebuild the bed, and a design variable can set it. `validate` reports a fixed dispersion velocity far from the adsorption velocity.
- **A process owns its flowsheet and cycle.** `Process(fs, cycle)` works on copies, so building it, or changing its design variables, leaves `fs` and `cycle` as they were, and several processes can be built from them.
- **Gradients with respect to design variables.** A variable's effect on each stage's settings and duration comes from central differences on the compiled cycle, which is cheap and, for the built-in steps, linear; the adjoints do the rest. Any numeric field of any step can then be a variable, including those of steps users define.
- **Energy counts compressors.** `performance` counts vacuum pumps on outlet devices below ambient pressure and compressors on inlet devices above it, as Haghpanah et al. do; the wet flue gas examples report the pumps alone too, as their paper does.
- **Roles are a fixed set**, unlike steps. The flowsheet diagram is printed text.
- **`simulate` extends `Jutul.simulate`**, so it does not clash with Jutul's.

Not done: plotting is still in `src/` (see the architecture proposal's next steps); files only hold isotherms given by component and the linear driving force; `history_matching.jl` still uses the old interface.
