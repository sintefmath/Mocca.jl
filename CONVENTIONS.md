# Conventions

The conventions Mocca follows in its inputs, names, signs and code, gathered
in one place so that the source reads quickly. Most are also stated in the
docstrings they belong to. [ARCHITECTURE.md](ARCHITECTURE.md) describes the
layers these sit in.

## Most likely to catch you out

1. **`AdsorbedLoading` is a copy of `AdsorbedConcentration`.** The primary
   variable is `AdsorbedConcentration`, the loading q. `AdsorbedLoading` is a
   secondary variable equal to it, which every time derivative uses so that
   sensitivities with respect to parameters work. `ColumnConservedEnergy`
   does the same for the temperature.
2. **Heats of adsorption differ in sign between isotherms.**
   `compute_enthalpy` returns a negative ΔH (adsorption is exothermic). The
   dual-site Langmuir `ΔUb` and `ΔUd` and the Toth–GAB `ΔH0` and `heats` are
   given negative. The Virial `ΔH` is given positive, as in the paper, and
   negated inside.
3. **A process works on copies.** `Process(fs, cycle)` deep-copies the
   flowsheet and the cycle, so changing them afterwards does not change the
   process. Change a process through its design variables. `p.flowsheet` and
   `p.cycle` are its copies.
4. **Solves start from the last steady state.** `cyclic_steady_state` starts
   from the last converged state and cycle Jacobian, kept in `p.warm`, so a
   result can depend on what was run before. Pass `warm_start = false` to
   start from the initial state.
5. **Pressures carry over from step to step.** A step without `from` starts
   where the step before left its bed. A variable on the evacuation
   pressure (`("evacuation", :to)`) also moves where the next pressurisation
   starts.
6. **Default boundaries follow the feed and the cycle.** A boundary set from
   a device role's defaults is updated whenever the feed or the cycle
   changes (`fs.default_boundaries`). One you set, with `boundary = ...` or
   `set_boundary!`, stays as it is.
7. **Heat capacities are given per kilogram** (J/(kg·K)), those of the gas
   included, but the energy balance works with molar ones: the state's
   `:C_pg` and `:C_pa` are in J/(mol·K), from the thermodynamics of the gas
   and of the adsorbed phase.

## Components and the key component

- **The order of the gas is the order of everything.** Component `i` is the
  `i`-th name in `Gas(names)`. Compositions `y`, loadings `q`, molar masses,
  isotherm vectors and the inventory closure all follow it.
- **The key component is the one to capture.** `key_component` (an index or a
  name, default 1) picks it, and by convention it comes first in the gas.
- **Performance is measured from the `key_*` integrals.** `key_extract` is the
  key component through the `:extract` devices, `total_extract` every
  component through them, and `key_feed` the key component through the
  `:feed` devices:
  - purity = `key_extract / total_extract`;
  - recovery = `key_extract / key_feed`;
  - specific energy is in kWh per tonne of key component in the extract;
  - productivity is in mol of key component per m³ of adsorbent per second.
- **The light product is the feed stripped of the key component**, cut to a
  trace of 1e-10 and the rest renormalised, unless the `Feed` is given one.
  Beds start full of it, and it is the composition at the boundary of
  `:product` and `:light_inlet` devices.
- **No mole fraction is ever zero.** `composition` gives the components it is
  not told about a trace of 1e-10. In the solver a mole fraction stays
  between `MIN_GAS_MOLEFRACTION` (1e-12) and one less that.
- **Material data are given by component name.** An `AbstractIsothermData`
  (`dual_site_langmuir`, `virial`, `toth_gab`, ...) is independent of the
  gas: components without parameters do not adsorb, and parameters for
  components the gas lacks are ignored. `build_isotherm` turns it into an
  `AbstractIsotherm` in the order of a particular gas.

## Units and signs

- **Plain SI units throughout:** Pa, K, m, s, mol, kg and J. The constants of
  `Mocca.Units` convert by multiplication, such as `0.1bar` (1e4 Pa) and
  `-36.6kJ_per_mol`, and `celsius(25)` is 298.15 K. Process files may give a
  number with its unit as a string, such as `"0.1 bar"`.
- **Docstrings give units in brackets** after a name, such as `to` [Pa].
- **Loadings are in mol per m³ of adsorbent** inside the solver. Isotherm
  data may be on a mass basis (mol/kg, `basis = :mass`), converted with the
  adsorbent density. Dual-site Langmuir data default to a volume basis, and
  Virial and Toth–GAB data to a mass basis.
- **Velocities:** `velocity` is interstitial; `superficial_velocity` and
  `flow_rate` (m³/s at the bed inlet) are the alternatives. The
  `dispersion_velocity` of a bed is interstitial.
- **The molar flow of a device is positive from its `:inlet` to its
  `:outlet`.** Feed and light inlet devices are connected to their bed by
  their outlet, so a positive flow enters the bed; product, extract and waste
  devices are connected by their inlet, so a positive flow leaves it. A link
  runs from its `from` port (its inlet) to its `to` port. Gas flowing
  backwards counts as negative in stream totals.
- **The inventory closure is positive when the beds gain** a component over
  the cycle.

## Ports, units and devices

- **A bed has ports `:bottom` (cell 1) and `:top` (its last cell).** In the
  VSA processes of the library, the feed enters at the bottom, the light
  product leaves at the top, and the extract is drawn by vacuum from the
  bottom, counter-current to the feed.
- **A device has ports `:inlet` (1) and `:outlet` (2)**, `DEVICE_PORTS`.
- **The sorbent of a moving bed moves down,** from its top cell to its bottom
  cell. `SorbentConnection(upstream, downstream)` passes the sorbent leaving
  the bottom of `upstream` into the top of `downstream`, through the ports
  `:sorbent_outlet` and `:sorbent_inlet`.
- **A port is a `PortRef(name, port)`**, or the tuple `(name, port)`, such as
  `(:Bed, :bottom)`.
- **Bed, unit with holdup and device are different things.** A bed is a unit
  that `holds_sorbent`; steps act on beds, and their adsorbent counts in the
  productivity. A unit with holdup (`holdup_units`) is any unit that is not a
  device. A device is a `FlowElement`: it has no holdup and owns the flow
  between its two ports.
- **The `Unit()` entity holds one value per unit**, such as the bed radius or
  the wall density, read in the equations as `state.Name[1]`. Mocca 0.1
  called it `Column`; `Column` is now the column geometry of the user
  interface.

## Names in a flowsheet

- **Units are CamelCase symbols:** `:Bed` for a single bed, `:A` and `:B` for
  two, `:SilicaGel` and `:Zeolite` for layers, `:Adsorption` and
  `:Desorption` for the sections of a moving bed.
- **Devices start with `V_`:** `:V_feed`, `:V_product`, `:V_vacuum` (extract),
  `:V_lpp` (light inlet, for light product pressurisation), `:V_eq`
  (equalisation link), `:V_link`, `:V_waste`, `:V_purge`. With several beds
  the bed follows, as in `:V_feed_A`.
- **A device's role decides what it does.** The roles (`DEVICE_ROLES`) are
  `:feed`, `:light_inlet`, `:product`, `:extract`, `:waste` and `:link`. The
  role sets:
  - which device a step opens;
  - which streams count as feed and as product in the performance;
  - the device's default boundary;
  - which work counts in the energy: vacuum pumps on outlet devices whose gas
    leaves below the ambient pressure, compressors on inlet devices whose gas
    enters above it, and blowers on links that raise the pressure.
- **A valve's default conductance** is that of the half cell at the end of
  the bed it is attached to (`default_conductance`), so that the valve adds no
  more resistance than the bed's end.

## Steps, stages and cycles

- **A step is written in process terms; a `Stage` is what it compiles to.**
  A step (`AbstractStep`) is one of `Pressurise`, `Adsorb`, `Blowdown`,
  `Evacuate`, `Purge`, `Equalise` and `Idle`, or one of your own. A stage is a
  duration and settings `device => (law = ..., inlet = ..., outlet = ...)`.
  Devices a stage does not name are closed. A unit with holdup may be named
  in a stage only to set its `jacket_temperature`.
- **The fields of the built-in steps:**
  - `to`: the pressure the step takes the bed to;
  - `from`: the pressure it starts from, `nothing` to carry on from the step
    before;
  - `through`: the end of the bed, `:top` or `:bottom`;
  - `source`: the gas a pressurisation uses, `:feed` or `:light_product`;
  - `bed`: the bed it acts on, `nothing` for the flowsheet's only bed;
  - `device`: a device to use instead of the one found by its role;
  - `ramp_rate` [1/s]: of the exponential pressure ramp, by default the
    process's (0.5).
- **Steps have default names:** `"pressurisation"` (or `"LPP"` from the light
  product), `"adsorption"`, `"blowdown"`, `"evacuation"`, `"purge"`,
  `"idle"` and `"equalisation A to B"`. With several beds the bed comes
  first, as in `"A evacuation"`. A `StepPlace` finds a step by this name or
  by its position.
- **A cycle of several sequences** runs them side by side, each labelled, such
  as one per bed. A shorter sequence is padded to the cycle time with its
  devices closed. A single sequence has the label `""`.
- **Time in a pressure profile** is time since the start of the stage.

## Processes

- **There are two kinds of `Process`.** A `CyclicProcess`, which
  `Process(fs, cycle)` builds, is operated by a cycle of steps. A
  `SteadyProcess`, such as a moving bed, keeps its device settings and runs
  to steady state.
- **Defaults of the options:** ambient pressure 101325 Pa and ambient
  temperature 298.15 K, pump efficiency 0.72, pressure ramp rate 0.5 per
  second, and the heat capacity ratio of the feed at its temperature
  (`pump_heat_capacity_ratio`); the presets and library processes set 1.4.
- **By default a bed starts** at the lowest pressure of its cycle and the
  feed temperature, full of light product, with the sorbent at equilibrium
  with the gas and the wall at the column's ambient temperature.
- **The compiled cycle is kept until it is needed again.** `_compiled(p)`
  compiles the cycle when it has to, `_refresh!(p)` after a change to the
  steps or the feed, and `_rebuild!(p)` sets up the model and parameters
  again after a unit was rebuilt.
- **Stage durations are differentiated through the `TimeScale` parameter**,
  which stretches the time steps of each stage. When a duration is a design
  variable, each stage keeps its number of time steps and they stretch with
  it (`TimeStepping(keep_step_pattern = ...)`); `reset_step_pattern!` takes
  the current durations as the new pattern.

## Performance and design variables

- **Performance is built from integrals, quantities and measures.**
  - Integrals are Jutul sum objectives over a cycle, or over one second at
    steady state, where they are rates: `key_extract`, `energy`, `heating`,
    ... (`performance_integrals`).
  - Quantities do not depend on the state: the `duration`, the adsorbent
    volume and mass, the sorbent inventory (`performance_quantities`).
  - Measures are functions `s -> ...` of both, where `s` holds the values of
    the integrals and quantities in one NamedTuple (`performance_measures`).
  - `add_integral!` and `add_measure!` add your own.
- **A place is a number in the description of a process, not in its state.**
  A place may change values only, never the structure of the model: the
  number of cells or stages and the isotherm are not places.
- **Places have tuple shorthands:**
  - `("evacuation", :to)`: a field of a step, by name or position;
  - `("A", "evacuation", :to)`: the same in the sequence labelled `"A"`;
  - `(:feed, :pressure)`: a field of the feed;
  - `(:parameter, :Bed, :AdsorbentDensity)`: a parameter of a unit;
  - `(:setting, :V_purge, :rate)` and `(:setting, :V_product, :outlet,
    :pressure)`: a device setting of a steady process;
  - `(:geometry, :Bed, :length)`: the geometry of a bed;
  - `ScaledPlace(place, factor)` sets a place to a multiple of the variable.
- **One variable may set several places**, tying them to one value. Read and
  set it as `p[:P_L]` and `p[:P_L] = 0.05bar`.
- **Variable names in the examples:** `:P_L` for the evacuation pressure,
  `:P_I` for the blowdown pressure, `:v_feed` for the feed velocity, and
  `:t_press`, `:t_ads`, `:t_blow` and `:t_evac` for the durations.
- **`optimise!` takes** `maximise = ...` or `minimise = ...`, targets as
  NamedTuples (`at_least = (purity = 0.9,)`, `at_most = ...`), a `scale` per
  measure, and `linear_constraints` as `(Dict(name => coefficient, ...),
  bound)`, meaning Σ coefficient·variable ≤ bound.

## Names in the Jutul model

- **The primary variables** of a bed are `:Pressure`, `:y` (mole
  fractions), `:AdsorbedConcentration`, `:Temperature` and, with a wall,
  `:WallTemperature`. A device has `:MolarFlow` and a copy of the state at
  each port: `:InletPressure`, `:InletTemperature`, `:InletComposition` and
  the same for `:Outlet`.
- **Per-component arrays are indexed `[component, cell]`.**
- **A few secondary variables have short names:** `:ΔH` (heat of
  adsorption), `:C_pa` and `:C_pg` (molar heat capacities of the adsorbed
  phase and the gas).
- **The equations** of a bed are `:mass_conservation`, `:mass_transfer`,
  `:energy_column` and `:energy_wall`, and of a device `:inlet_state`,
  `:outlet_state` and `:device_law`.
- **Forces:** a unit with holdup takes `jacket_temperature` and `time_scale`,
  a device `law`, `inlet` and `outlet`.
- **Parameters are CamelCase,** and a process overrides them per unit:
  `parameters = Dict(:Bed => Dict(:AdsorbentDensity => 1200.0))`.

## Convergence

- **The change over a cycle is tested with mixed absolute and relative
  tolerances** (`cycle_change`). Each variable is measured against its own
  scale, per component for mole fractions and loadings (at least 1e-6, and
  0.1% of the saturation capacity), less what the time steps resolve at the
  nonlinear tolerance.
- **Default tolerances:** Newton's method converges the cycle to 1e-9 with
  time steps converged to 1e-8; Anderson acceleration to 1e-6 with time
  steps converged to 1e-6. A device's equations converge 1000 times tighter
  than the rest.
- **The inventory closure (`closure`) is reported, but does not decide
  convergence.**

## Code

- **Each layer uses only the layers below it:** Framework, Physics,
  Balances, then the units (Devices, FixedBeds, MovingBeds), Networks,
  Operation and Interface. Legacy and Plotting sit beside them.
  - A module brings in a layer whole, `using ..Balances`, and imports only
    helpers whose names start with `_` by name.
  - A method added to a function of another module is written with the
    function's full name, `Framework.ports(model::FlowDeviceModel) = ...`.
  - Each source folder is named after its module, in lower case.
  - Unit modules are plural (`FixedBeds`), so that `Mocca.FixedBed` stays
    the type.
- **A name starting with `_` is internal:** not exported, and free to change.
- **Errors in the input** are thrown as `input_error(location, field,
  message)`, a `MoccaInputError` that says where the problem is.
  `_at(location) do ... end` turns an error from a lower layer into one.
  A `SimulationFailure` means the solver gave up, not that the input or the
  code is wrong.
- **Name patterns:**
  - `setup_*` builds Jutul objects, and `compute_*` evaluates physics;
  - a function ending in `!` changes its first argument;
  - `*Objective` types are Jutul objectives, called as
    `(model, state, dt, step_info, forces)`;
  - `*Spec` records how a unit was built, to rebuild it and write it to a
    file;
  - `*Data` holds material data by component name.
- **Docstrings** start with the call signature and `→` the result, and give
  units in brackets. Comments say why the code does what it does.
- **Library files** are named `<material>_<author>_<year>` for adsorbents,
  such as `zeolite_13x_haghpanah_2013`, and `<author>_<year>_<process>` for
  processes, such as `haghpanah_2013_four_stage_vsa`.
