```@meta
CurrentModule = Mocca
```

# Model equations

A fixed bed in Mocca is a one-dimensional, axially dispersed plug flow model
of a packed column, as in Haghpanah et al. (2013). The gas is ideal, flow
follows Darcy's law, mass transfer to the sorbent is a linear driving force,
and the energy balance covers the gas, the sorbent and the column wall. The
beds are connected to the feed, products and each other through flow devices
(valves and flow controllers), whose settings change from step to step of the
cycle.

All quantities are in SI units. The symbols are listed in the
[nomenclature](#Nomenclature) at the end of this page, with the inputs they
come from.

## Fixed bed

### Component mass balance

For each component ``i`` of the gas,

```math
\varepsilon \frac{\partial c_i}{\partial t}
+ \frac{\partial (u\, c_i)}{\partial z}
- \frac{\partial}{\partial z}\left(\varepsilon D_L C \frac{\partial y_i}{\partial z}\right)
+ (1-\varepsilon)\frac{\partial q_i}{\partial t} = 0,
\qquad c_i = y_i C = \frac{y_i P}{RT}.
```

``q_i`` is the loading in mol per m³ of adsorbent (particles), so
``(1-\varepsilon)`` converts it to bed volume.

### Momentum: Darcy's law

The superficial velocity follows Darcy's law with the Kozeny–Carman
(Blake–Kozeny) permeability of a packed bed,

```math
u = -\frac{k}{\mu}\frac{\partial P}{\partial z},
\qquad
k = \frac{\varepsilon^3 d_p^2}{150\,(1-\varepsilon)^2}.
```

This is the viscous (laminar) term of the Ergun equation; the inertial term is
not included. The interstitial velocity is ``v = u/\varepsilon``. Give
`permeability` to [`add_bed!`](@ref) to override ``k``.

### Axial dispersion

```math
D_L = 0.7\, D_m + 0.5\, v_0\, d_p,
```

where ``D_m`` is the molecular diffusivity of the gas (`diffusivity` of
[`Gas`](@ref), that of the first component if given per component) and ``v_0`` the interstitial feed velocity of the bed's
adsorption step. ``D_L`` is evaluated once, when the process is built, and is
constant along the bed and through the cycle. Set `dispersion_velocity` or
`dispersion` in [`add_bed!`](@ref) to fix it yourself.

### Adsorption rate: linear driving force

```math
\frac{\partial q_i}{\partial t} = k_i \left(q_i^* - q_i\right),
\qquad
k_i = \frac{c_i}{q_i^*}\,\frac{15\,\varepsilon_p D_{m,i}}{\tau\, r_p^2},
\qquad r_p = \frac{d_p}{2}.
```

This is macropore diffusion control: the particle porosity ``\varepsilon_p``,
tortuosity ``\tau`` and diameter ``d_p`` of the [`Adsorbent`](@ref) set the rate.

### Isotherm: dual-site Langmuir

```math
q_i^* = \frac{q_{sb,i}\, b_i c_i}{1 + \sum_j b_j c_j}
      + \frac{q_{sd,i}\, d_i c_i}{1 + \sum_j d_j c_j},
\qquad
b_i = b_{0,i}\, e^{-\Delta U_{b,i}/RT},
\quad
d_i = d_{0,i}\, e^{-\Delta U_{d,i}/RT}.
```

The affinities are on a concentration basis (m³/mol), and the saturation
capacities in mol per m³ of adsorbent, or in mol/kg with `basis = :mass`,
which Mocca multiplies by the adsorbent density. The single-site
[`langmuir`](@ref) isotherm has ``q_{sd,i} = 0``. Check an isotherm against its
source with [`equilibrium_loading`](@ref) and [`isotherm_table`](@ref), which
report mol/kg.

The heat of adsorption of each component is constant, the average of its two
site energies weighted by their capacities at the reference temperature
``T_0`` (`reference_temperature`, default 298.15 K):

```math
\Delta H_i = \frac{q_{sb,i}\left(\Delta U_{b,i} - RT_0\right) + q_{sd,i}\left(\Delta U_{d,i} - RT_0\right)}{q_{sb,i} + q_{sd,i}}.
```

### Energy balance of the column

Written in terms of temperature and pressure:

```math
\begin{aligned}
&\left[(1-\varepsilon)\left(\rho_s c_{ps} + \tilde c_{pa} \sum_i q_i\right)\right]\frac{\partial T}{\partial t}
+ \varepsilon\,\frac{\tilde c_{pg}}{R}\frac{\partial P}{\partial t}
+ \frac{\tilde c_{pg}}{R}\frac{\partial (uP)}{\partial z}
- \frac{\partial}{\partial z}\left(K_z \frac{\partial T}{\partial z}\right) \\
&\qquad + (1-\varepsilon)\sum_i \left(\tilde c_{pa} T + \Delta H_i\right)\frac{\partial q_i}{\partial t}
+ \frac{2 h_{in}}{r_{in}}\left(T - T_w\right) = 0,
\end{aligned}
```

with molar heat capacities of the gas and the adsorbed phase from the
per-mass values of each component,

```math
\tilde c_{pg} = \bar M \sum_i y_i\, c_{pg,i},
\qquad
\tilde c_{pa} = \bar M \sum_i y_i\, c_{pa,i},
\qquad
\bar M = \sum_i y_i M_i.
```

``\Delta H_i < 0``, so adsorption heats the bed.

### Wall energy balance

```math
\rho_w c_{pw}\frac{\partial T_w}{\partial t}
= K_w \frac{\partial^2 T_w}{\partial z^2}
+ \frac{2 r_{in} h_{in}}{r_{out}^2 - r_{in}^2}\left(T - T_w\right)
- \frac{2 r_{out} h_{out}}{r_{out}^2 - r_{in}^2}\left(T_w - T_a\right).
```

At both ends the wall conducts heat to the surroundings at ``T_a``.

### Thermal options

The `thermal` keyword of [`add_bed!`](@ref) chooses the energy model:

| `thermal` | Energy model |
|---|---|
| `WithWall()` (default) | Column and wall energy balances, as above |
| `Adiabatic()` | Column energy balance without the wall term (``h_{in} = 0``) |
| `Isothermal(T)` | Temperature fixed at `T`, no energy balance |

## Flow devices and boundary conditions

Each end of a bed (`:bottom` at ``z = 0``, `:top` at ``z = L``) connects to one
or more flow devices. A device has an inlet and an outlet port and carries a
molar flow ``F`` [mol/s], positive from inlet to outlet. The port attached to
the bed takes the pressure, temperature and composition of the end cell. The
other port is either another bed (a link) or a boundary at a prescribed
pressure, temperature and composition. The flow enters or leaves the bed's
end cell carrying the composition and enthalpy ``F\,\tilde c_{pg} T`` of the
upstream port.

In each step of the cycle every device is closed, a valve or a flow
controller:

| Law | Flow |
|---|---|
| Closed | ``F = 0`` |
| [`LinearValve`](@ref) | ``F = C_v\,(P_{in} - P_{out})\,\dfrac{P_{up}}{R\,T_{up}}`` |
| [`VolumetricFlow`](@ref) (flow controller) | ``F = Q\,\dfrac{P_{out}}{R\,T_{in}}`` |

Here ``C_v`` is the valve conductance [m³/(s·Pa)] and "up" is the upstream port.
A feed given as an interstitial velocity ``v`` is a flow controller with
``Q = \varepsilon v A``.

By default a valve has the conductance of half a cell at the bed end it is
attached to, ``C_v = kA/(\mu\,\Delta z/2)``, so that it adds no more resistance
than the bed itself and the bed end follows the pressure on the far side of
the valve. Give `conductance` to a step to model a smaller valve.

During pressurisation, blowdown and evacuation, the pressure on the far side
of the valve moves exponentially from where the step starts to where it ends:

```math
P(t) = P_{end} + \left(P_{start} - P_{end}\right) e^{-\lambda t},
```

with ``\lambda`` the `ramp_rate` (default 0.5 s⁻¹) of the step or the
[`Process`](@ref).

For the four-step VSA cycle this gives the boundary conditions of Haghpanah et
al. (2013):

| Step | Bottom (``z = 0``) | Top (``z = L``) |
|---|---|---|
| Pressurisation | feed valve, ``P`` ramps ``P_L \to P_H`` | closed |
| Adsorption | feed at velocity ``v_{feed}`` | product valve at ``P_H`` |
| Blowdown | closed | product valve, ``P`` ramps ``P_H \to P_I`` |
| Evacuation | vacuum valve, ``P`` ramps ``P_I \to P_L`` | closed |

## Performance

[`performance`](@ref) measures a process over its last simulated cycle, from
the streams through its devices. The device roles (see [`add_device!`](@ref))
say which streams count. For the key component ``k`` (the first component of
the gas unless `key_component` is given to [`Process`](@ref)), with ``N_i^{r}``
the moles of component ``i`` that passed through all devices of role ``r``
over the cycle:

```math
\text{purity} = \frac{N_k^{extract}}{\sum_i N_i^{extract}},
\qquad
\text{recovery} = \frac{N_k^{extract}}{N_k^{feed}},
\qquad
\text{productivity} = \frac{N_k^{extract}}{V_{ads}\; t_{cycle}},
```

with ``V_{ads}`` the volume of adsorbent, ``(1-\varepsilon)`` times the bed
volume, summed over the beds. Productivity is in mol/(m³ adsorbent · s);
multiply by ``M_k \cdot 86400/1000`` for tonnes per m³ per day.

The energy counts a vacuum pump on every outlet device (roles `:product`,
`:extract`, `:waste`) while it discharges below ambient pressure, and a
compressor on every inlet device (`:feed`, `:light_inlet`) while it delivers
above ambient pressure. Both compress adiabatically between the ambient pressure
``P_{amb}`` and the device's pressure ``P``:

```math
E = \int \frac{F R T}{\eta}\,\frac{\gamma}{\gamma - 1}
\left[\left(\frac{P_{high}}{P_{low}}\right)^{(\gamma-1)/\gamma} - 1\right] dt,
\qquad
\text{specific energy} = \frac{E}{N_k^{extract} M_k},
```

reported in kWh per tonne of key component. ``P_{amb}``, ``\eta`` and
``\gamma`` are the `ambient_pressure` (101325 Pa), `pump_efficiency` (0.72) and
`heat_capacity_ratio` (1.4) of the [`Process`](@ref). Haghpanah et al. count
the vacuum pumps only; [`vacuum_pump_energy`](@ref) gives the energy of one
device, to report it their way.

## Cyclic steady state

The cycle maps the state of the beds at its start, ``x``, to the state at its
end, ``\Phi(x)``. [`cyclic_steady_state`](@ref) solves ``x = \Phi(x)``:

| `method` | How | Converged when |
|---|---|---|
| `:newton` (default) | Anderson acceleration, then Newton's method on ``x - \Phi(x)`` | change per cycle below `tol` = 10⁻⁹ |
| `:anderson` | Repeated cycles with Anderson acceleration | change per cycle below `tol` = 10⁻⁶ |
| `:mass_balance` | Repeated cycles | key component in and out agree within 0.5% for 5 cycles in a row, after at least 50 (Haghpanah et al. 2013; Krishnamurthy et al. 2014) |

The change per cycle is the largest relative change of any variable of any
bed (see [`cycle_change`](@ref)). `:newton` and `:anderson` require every
variable to settle, including slowly accumulating components such as water; a
mass balance on one component does not.

## Numerical method

Each bed is divided into `ncells` equal cells (30 by default). Fluxes between
cells are two-point, with upwinded advection, and time stepping is fully
implicit (backward Euler), solved by Newton's method with exact derivatives
from automatic differentiation. Each step of the cycle starts with time steps
of `first_dt` that double up to `max_dt` (see [`TimeStepping`](@ref)), which
carries the solver through the abrupt start of a step. Gradients for
optimisation come from the adjoint of the same discrete equations.

## Nomenclature

| Symbol | Meaning | Unit | Set by |
|---|---|---|---|
| ``A`` | Column cross-section, ``\pi r_{in}^2`` | m² | [`Column`](@ref) `radius` |
| ``b_0, d_0`` | Affinity at infinite temperature, sites b and d | m³/mol | isotherm `b0`, `d0` |
| ``c_i,\ C`` | Gas concentration of component ``i``, total | mol/m³ | state |
| ``c_{pa,i}`` | Heat capacity of adsorbed component ``i`` | J/(kg·K) | [`Adsorbent`](@ref) `adsorbed_heat_capacity` |
| ``c_{pg,i}`` | Heat capacity of gas component ``i`` | J/(kg·K) | [`Gas`](@ref) `heat_capacities` |
| ``c_{ps}`` | Heat capacity of the adsorbent | J/(kg·K) | [`Adsorbent`](@ref) `heat_capacity` |
| ``c_{pw}`` | Heat capacity of the wall | J/(kg·K) | [`Wall`](@ref) `heat_capacity` |
| ``C_v`` | Valve conductance | m³/(s·Pa) | step `conductance` |
| ``d_p`` | Particle diameter | m | [`Adsorbent`](@ref) `particle_diameter` |
| ``D_L`` | Axial dispersion coefficient | m²/s | computed, or [`add_bed!`](@ref) `dispersion` |
| ``D_m`` | Molecular diffusivity | m²/s | [`Gas`](@ref) `diffusivity` |
| ``F`` | Molar flow through a device | mol/s | state |
| ``h_{in}, h_{out}`` | Heat transfer coefficient, bed to wall and wall to ambient | W/(m²·K) | [`Column`](@ref) `inner_heat_transfer`, `outer_heat_transfer` |
| ``k`` | Bed permeability | m² | computed, or [`add_bed!`](@ref) `permeability` |
| ``K_z`` | Axial thermal conductivity of the bed | W/(m·K) | [`Column`](@ref) `axial_conductivity` |
| ``K_w`` | Thermal conductivity of the wall | W/(m·K) | [`Wall`](@ref) `conductivity` |
| ``L`` | Column length | m | [`Column`](@ref) `length` |
| ``M_i`` | Molar mass | kg/mol | [`Gas`](@ref) `molar_masses` |
| ``P`` | Pressure | Pa | state |
| ``q_i,\ q_i^*`` | Loading and equilibrium loading | mol/m³ adsorbent | state |
| ``q_{sb}, q_{sd}`` | Saturation capacity, sites b and d | mol/m³ adsorbent | isotherm `qsb`, `qsd` |
| ``r_{in}, r_{out}`` | Inner and outer radius of the column wall | m | [`Column`](@ref) `radius`, `outer_radius` |
| ``R`` | Gas constant, 8.314 | J/(mol·K) | |
| ``T,\ T_w,\ T_a`` | Gas and sorbent, wall and ambient temperature | K | state; [`Column`](@ref) `ambient_temperature` |
| ``u,\ v`` | Superficial and interstitial velocity | m/s | state; step `velocity` |
| ``y_i`` | Gas mole fraction | – | state |
| ``\Delta U_b, \Delta U_d`` | Internal energy of adsorption, sites b and d | J/mol | isotherm `ΔUb`, `ΔUd` |
| ``\varepsilon`` | Bed porosity (void fraction between particles) | – | [`Column`](@ref) `porosity` |
| ``\varepsilon_p`` | Particle porosity | – | [`Adsorbent`](@ref) `particle_porosity` |
| ``\gamma,\ \eta`` | Heat capacity ratio and efficiency of pumps and compressors | – | [`Process`](@ref) |
| ``\lambda`` | Pressure ramp rate | 1/s | step or [`Process`](@ref) `ramp_rate` |
| ``\mu`` | Gas viscosity | Pa·s | [`Gas`](@ref) `viscosity` |
| ``\rho_s,\ \rho_w`` | Density of the adsorbent and of the wall | kg/m³ | [`Adsorbent`](@ref) `density`, [`Wall`](@ref) `density` |
| ``\tau`` | Tortuosity | – | [`Adsorbent`](@ref) `tortuosity` |
