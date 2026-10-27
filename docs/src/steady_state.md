```@meta
CurrentModule = Mocca
```

# Cyclic steady state

A cyclic process never settles to a constant state: it repeats a cycle of
steps. After enough cycles the state of the beds at the end of a cycle equals
the state at its start. That is the cyclic steady state, and performance is
measured over that cycle. [`cyclic_steady_state`](@ref) finds it:

```julia
result = cyclic_steady_state(process)
result.converged, result.iterations
```

The cycle maps the state of the beds at its start, ``x``, to the state at its
end, ``\Phi(x)``, and the steady state solves ``x = \Phi(x)``. The change over
a cycle measures how far ``\Phi(x)`` is from ``x``, as described under
[How convergence is judged](@ref) below. The page
[How a VSA cycle works](@ref) shows the approach to the steady state of one
process.

## Methods

| `method` | How | Converged when |
|---|---|---|
| `:newton` (default) | A few plain cycles, then Newton's method on ``x - \Phi(x)``, every cycle with steps converged to 10⁻⁸ | change per cycle below `tol` = 10⁻⁹ |
| `:anderson` | Repeated cycles with Anderson acceleration, steps converged to 10⁻⁶ | change per cycle below `tol` = 10⁻⁶ |
| `:mass_balance` | Repeated cycles | key component in and out agree within 0.5% for 5 cycles in a row, after at least 50 (Haghpanah et al. 2013; Krishnamurthy et al. 2014) |

`:newton` and `:anderson` require every variable to settle, including slowly
accumulating components such as water. A mass balance on one component does
not, and can stop while the water is still far from steady.

## How convergence is judged

Each variable of each bed, at each cell, has settled when its change over a
cycle is within the tolerance relative to its scale plus what the time steps
resolve, as the mixed absolute and relative tolerances of an ODE solver:

```math
|\Phi(x)_i - x_i| \le \mathrm{tol} \cdot s_i + a_i
```

The scale ``s_i`` is the largest magnitude of the variable over the bed at the
start of the cycle: per component for mole fractions and loadings, so that a
minor component is measured against its own level, and over the whole bed for
pressures and temperatures. The resolution ``a_i`` comes from the nonlinear
tolerance of the steps, `nonlinear_tolerance`. Each step converges the mass
balance of a cell to that fraction of the cell's holdup, and its energy
balance to that error in kelvin, so a step resolves

| Variable | Resolution ``a_i`` |
|---|---|
| mole fraction | the nonlinear tolerance |
| loading | the nonlinear tolerance times the cell's total loading |
| pressure | the nonlinear tolerance times the pressure |
| temperature | the tolerance of its energy balance, in K: the nonlinear tolerance, but at least 10⁻⁷ K for the gas |

A change within that is the noise of the steps, and asking for less would ask
for more than the simulation gives. The change reported for a cycle, in the
progress and in the result's `history`, is the largest
``\max(|\Phi(x)_i - x_i| - a_i, 0) / s_i``: the cycle passes when it is below
`tol`, and it is zero once nothing changes by more than the steps resolve.
Far from the steady state it is the largest relative change.

Every component's mole fraction is checked. The simulation solves for all
but the last, which follows from the others, so the last changes by minus
the sum of their changes; with the resolution in the test, that does not
hold Newton's method back when the last component is a trace, and the
result does not depend on the order of the components. Newton's method
weighs its residual the same way, so that it reduces what the test measures.

All the cycles of a solve are simulated at the same nonlinear tolerance, so
that each phase works with the same cycle map: with `:newton`, the plain
cycles and any Anderson acceleration before Newton's method run at Newton's
10⁻⁸, which is also that of the adjoints for gradients. Pass
`nonlinear_tolerance` to change it, for example to converge more tightly.

The result also has the inventory closure of the last cycle, `closure`: for
each component, the change of the moles the beds hold, in the gas and on the
sorbent, over the cycle, as a fraction of the moles of it fed through the
inlet devices. It is printed at the end of the progress, here for the wet
flue gas LPP VSA:

```
Inventory change over the last cycle, per mole fed: CO2 -8.68e-10, N2 1.59e-16, H2O 3.89e-10
```

It does not depend on the order of the components either, and it shows a
component still accumulating, but it does not decide convergence.

## How Newton's method gets there

Three plain cycles (`handover`), each from the end of the last, let the
fastest changes die out. Newton's method then solves ``x - \Phi(x) = 0``. Its
Jacobian, ``I - \partial\Phi/\partial x``, comes from forward sensitivities
over one cycle, every column in one sweep, so a Jacobian costs about as much
as a cycle on a coarse grid.

Far from the steady state a Newton step can make the change larger before the
next ones converge. A step that makes it more than 100 times larger, or from
which the cycle fails to simulate, is cut back, up to six times. If Newton's
method still does not converge, Anderson acceleration takes the change to
10⁻¹ and then 10⁻³, and Newton's method starts again from each. With
`handover = 1e-3`, Anderson acceleration takes the change to `1e-3` before
Newton's method starts, as Mocca 0.2.0 did.

With `warm_start` (the default), a new run starts from the last steady state
found for the process, and its Jacobian. After a small change of a setting it
converges in a few iterations, which is what makes sweeps and optimisation
quick.

## Fine grids

Two more things keep Newton's method fast on fine grids:

| Option | Default | Effect |
|---|---|---|
| `krylov` | `:auto` | Solves each Newton step by GMRES, one forward sensitivity per product, instead of building the full Jacobian. `:auto` does so when the beds have more than 800 variables (a bed of about 130 cells with two components, 100 with three); `true` or `false` fixes the choice, and a number moves the threshold. |
| `coarse_cells` | 40 | On a bed of more cells, finds the steady state first on beds of half as many cells, down to 40 or fewer, converged to 10⁻⁶, and starts each finer grid from the coarser one's, interpolated. `nothing` starts every grid from its initial state. |

The coarser grids only give the start; the steady state is that of the
process as set up. If Newton's method does not converge from a coarser
grid's steady state, the fine grid is started from its initial state.

## Options

| Keyword | Default | Meaning |
|---|---|---|
| `method` | `:newton` | See Methods above |
| `tol` | 10⁻⁹ (`:newton`), 10⁻⁶ (`:anderson`) | Change over a cycle to stop at |
| `max_cycles` | 500 | Most cycles to run |
| `warm_start` | `true` | Start from the last steady state of the process |
| `handover` | 3 | Plain cycles before Newton's method; a value below 1 is a change for Anderson acceleration to reach instead |
| `krylov` | `:auto` | GMRES or the full Jacobian for each Newton step |
| `coarse_cells` | 40 | Grid sequencing above this many cells; `nothing` turns it off |
| `nonlinear_tolerance` | 10⁻⁸ (`:newton`), 10⁻⁶ (`:anderson`) | Tolerance of the steps of every cycle, and so the resolution in the test |
| `progress` | on in the REPL | Print each cycle as it runs; an `IO` to print to |

The full list, and the lower-level solvers, are in the docstrings of
[`cyclic_steady_state`](@ref), [`newton_cyclic_steady_state`](@ref),
[`simulate_to_cyclic_steady_state`](@ref) and
[`simulate_until_mass_balance`](@ref).

## Following the progress

In the REPL, [`cyclic_steady_state`](@ref) prints each cycle as it runs, with
a line and a bar for each phase. In a terminal the lines are rewritten in
place. While Newton's method runs:

```
Cyclic steady state of four-stage VSA: 3 plain cycles, then Newton's method to 1.0e-09
  Cycles          3 cycles      change 6.51e-01            ▕████████████▏ 100%  done, on to Newton  0.5 s
  Newton        iteration   5   change 6.31e-02 → 1.0e-09  ▕█           ▏   8%  1.7 s
```

and when it is done:

```
Cyclic steady state of four-stage VSA: 3 plain cycles, then Newton's method to 1.0e-09
  Cycles          3 cycles      change 6.51e-01            ▕████████████▏ 100%  done, on to Newton  0.5 s
  Newton         13 iterations  change 1.96e-10 < 1.0e-09  ▕████████████▏ 100%  converged  4 Jacobians built  3.2 s
Converged after 17 cycles simulated in 4.1 s
```

A step that was cut back is counted on the line. With GMRES, the line shows
how far the solve of the current step has come:

```
  Newton        iteration   0   change 3.19e-01 → 1.0e-09  ▕            ▏   0%  GMRES: product 12, residual 3.4e-02 → 1.0e-02
```

With grid sequencing, each grid has its own lines (here with
`krylov = true`, at 160 cells):

```
Cyclic steady state of four-stage VSA: 3 plain cycles, then Newton's method on 40 cells to 1.0e-06; from there Newton's method on 80 cells to 1.0e-06 and on 160 cells to 1.0e-09
  40 cells: 3 plain cycles, then Newton's method
  Cycles          3 cycles      change 4.49e-01            ▕████████████▏ 100%  done, on to Newton  1.1 s
  Newton         10 iterations  change 1.18e-08 < 1.0e-06  ▕████████████▏ 100%  converged  273 GMRES products  8.4 s
  80 cells: Newton's method from the steady state on 40 cells
  Newton          2 iterations  change 5.89e-07 < 1.0e-06  ▕████████████▏ 100%  converged  76 GMRES products  3.0 s
  160 cells: Newton's method from the steady state on 80 cells
  Newton          3 iterations  change 7.55e-11 < 1.0e-09  ▕████████████▏ 100%  converged  114 GMRES products  7.1 s
Converged after 21 cycles simulated in 20 s
```

Each bar shows how far the smallest change so far has come from the first
cycle of its phase to the phase's tolerance, on a log scale, and the estimate
of what is left is for that phase. From the last steady state of the process,
Newton's method starts at once, with a single line. With
`method = :mass_balance` the line shows the balance error of the key
component instead, and how many cycles in a row it has been below the
tolerance.

Pass `progress = false` to turn it off, or `progress = true` to have it in a
script, as the examples do. Outside a terminal, such as in a log file, each
cycle gets its own line. It is off when [`optimise!`](@ref) finds the steady
state of each design.
