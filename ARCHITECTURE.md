# Mocca 0.2 architecture

Mocca 0.2 replaces the single, hard-wired adsorption column of 0.1 with a
flowsheet of units joined by flow devices, built on
[Jutul](https://github.com/sintefmath/Jutul.jl). The same model runs a one-bed
VSA cycle or a multi-bed process. It finds the cyclic steady state and gives
adjoint gradients of the performance with respect to operating settings, stage
durations and model parameters.

## Layers

Each layer is a submodule of `Mocca` and only uses the layers below it. Most
users work in the top layer, and `using Mocca` exports what they need. The
layers underneath are open for new units, steps and physics.

| Layer | Module | Source | Role | Main types |
|---|---|---|---|---|
| **Interface** | `Mocca.Interface` | `src/interface` | Describe a process in process terms, in SI units. Every part can be saved to and read from JSON, and a library holds the published processes. | `Gas`, `Feed`, `Adsorbent`, `Column`, `Flowsheet`, `add_bed!`, `add_device!`, `Cycle`, `Pressurise`, `Adsorb`, `Evacuate`, `Equalise`, `Process`, `add_variable!`, `optimise!` |
| **Process** | `Mocca.Operation` | `src/operation` | Stages and time steps, cyclic steady state, performance metrics, and adjoint gradients, including at cyclic steady state. | `Stage`, `setup_schedule`, `parallel_stages`, `newton_cyclic_steady_state`, `force_gradients`, `purity`, `recovery` |
| **Flowsheet** | `Mocca.Networks` | `src/networks` | A network of named units, the connections between their ports and the boundary conditions, turned into one Jutul `MultiModel`. Connections become cross terms. | `UnitNetwork`, `connect!`, `StreamConnection`, `PortStateCT`, `StreamCT` |
| **Units** | `Mocca.FixedBeds`, `Mocca.Devices` | `src/fixedbeds`, `src/devices` | One module per kind of unit. A fixed bed is a column of cells with holdup, along its axis and optionally in rings across its radius. Flow devices are one-cell models that own the molar flow between two ports. | `FixedBed`, `FlowDevice`, `LinearValve`, `VolumetricFlow`, `PortCondition`, `ExponentialRamp` |
| **Balances** | `Mocca.Balances` | `src/balances` | Shared material and energy balances: gas mass, adsorbed loading, gas energy, wall energy. Each unit adds its physics as a tuple of source terms. Scale-free convergence checks, and the setup every bed of sorbent shares. | `SorptionUptake`, `HeatOfSorption`, `WallExchange`, `setup_process_state` |
| **Physics** | `Mocca.Physics` | `src/physics` | Pluggable material models behind small interfaces. | `AbstractIsotherm`, `DualSiteLangmuir`, `AbstractMassTransfer`, `LinearDrivingForce` |
| **Core** | `Mocca.Framework` | `src/framework` | The system hierarchy, the core variables, thermal models, and the functions every unit implements. | `MoccaSystem`, `SorbentBed`, `ports`, `stream_equations`, `WithWall`, `Adiabatic`, `Isothermal` |
| **Jutul** | dependency | | Discretisation, automatic differentiation, Newton solver and adjoints. | |

Beside the layers, `Mocca.Units` (`src/Units.jl`) holds the unit conversions,
`Mocca.Legacy` (`src/legacy`) reads Mocca 0.1 inputs and `Mocca.Plotting`
(`src/plotting`) draws the results. Each folder is named after its module.

A module brings in the layers it builds on whole, with `using ..Balances`, and
exports what the layers above use. Only helpers whose names start with `_` are
imported by name. A method added to a function of another module is written
with the function's full name, such as `Framework.ports(model::FlowDeviceModel)`.
Without it, Julia either rejects the method or quietly defines a new function
of the same name, whose methods are never called.

## From a cycle to results

1. The `Cycle` of steps is compiled into `Stage`s. Each step opens the devices
   it needs, found by their role (feed, product, extract, link), and sets their
   flow law and port conditions. Pressure ramps carry on from where the step
   before left the bed.
2. The stages become Jutul forces and time steps. The flowsheet becomes a
   `MultiModel` with one model per bed and per device.
3. The cyclic steady state is found by repeating cycles with Anderson
   acceleration, then Newton's method on `x = Φ(x)` with an adjoint Jacobian.
4. Purity, recovery, productivity and energy are measured from the streams
   through the devices.
5. Adjoint gradients of these metrics with respect to named design variables
   drive the optimiser.

## Key design choices

- **Devices own the flow, units own the holdup.** A bed never sets its own
  boundary conditions; it connects through devices, as wells couple to a
  reservoir in JutulDarcy. Each component and the stream enthalpy are conserved
  across every connection.
- **A port is a layer of cells.** A unit meshed across its cross-section,
  such as a bed in rings, connects through the same devices as a column: the
  stream is shared between the cells of the port by their cross-section, and
  each gives up its own gas, mixed in the device.
- **Boundary conditions are device settings.** The fixed set of stage
  boundary-condition types is gone. A stage is a valve law plus the pressure,
  temperature and composition at a device's open port.
- **One set of conservation laws, physics as source terms.** A new unit lists
  its sorption, heat and wall terms rather than rewriting the residuals. The
  term tuples are known at compile time, so unused terms cost nothing.
- **Everything is a differentiable parameter.** Device settings, stage
  durations and unit parameters such as axial dispersion all get adjoint
  gradients, also at cyclic steady state.
- **Convergence is scale-free.** Mass errors are measured relative to the
  cell's holdup and energy errors as a temperature change in K. Cycle-to-cycle
  changes are measured per component, so trace components count.

## Extension points

| To add | Implement |
|---|---|
| An isotherm | `compute_equilibrium`, `compute_enthalpy`, and optionally `saturation_capacity` |
| A mass transfer model | `compute_mass_transfer_rate` |
| A source term | `mass_source`, `energy_source`, `wall_source` or `heat_capacity`, then list it for the unit |
| A unit | A `MoccaSystem` with `ports`, `stream_equations` and `holds_sorbent`, in a module of its own like `Mocca.FixedBeds` |
| A device law | `device_flow` for a subtype of `AbstractDeviceLaw` |
| A step | `step_settings` for a subtype of `AbstractStep`, and `register_step!` for files |
| A connection | `connection_ports`, `check_connection` and `add_cross_terms!` |

## Changes from 0.1

| Mocca 0.1 | Mocca 0.2 |
|---|---|
| One column with its constants struct | Beds, devices and links in a flowsheet |
| A boundary-condition type per stage | Steps compiled to device settings |
| Stage times as the only cycle input | Cycles of steps, including concurrent multi-bed steps |
| Heat of sorption of the first component only | Every component counts |
| Cycles repeated until converged | Anderson acceleration, then Newton's method |
| Gradients of a single simulation | Gradients at cyclic steady state, including stage durations |

See [NEWS.md](NEWS.md) for the full release notes.

## In use

```julia
using Mocca
using Mocca.Units

# Haghpanah et al. (2013) four-stage VSA from the library
process = read_process("haghpanah_2013_four_stage_vsa")
process[:P_L] = 0.05bar                 # a design variable
result = cyclic_steady_state(process)
performance(result)                     # purity, recovery, kWh/t, mol/(m³·s)
performance_gradient(process)           # adjoint gradients of each metric
```
