# Utility functions for mass and energy balance checks

# --- Private helpers ---

"""
Compute BC transmissibility at a boundary cell directly from model domain data.
"""
function _bc_trans_from_model(model::AdsorptionModel, cell)
    k = model.data_domain[:permeability][cell]
    r_in = first(model.data_domain[:r_in, Column()])
    dx = model.data_domain[:dx][cell]
    A = π * r_in^2
    return k * A / (dx / 2.0)
end

"""
Compute the column face area from the mesh geometry.
"""
function _column_face_area_from_model(model::AdsorptionModel)
    g = Jutul.physical_representation(model.data_domain)::Jutul.CartesianMesh
    return g.deltas[2] * g.deltas[3]
end

"""
Return the force at step `i` from a vector of forces or a single force applied uniformly.
"""
_get_force_at_step(forces::AbstractVector, i) = forces[i]
_get_force_at_step(forces, i) = forces

"""
Compute total moles of component `comp` in the column (fluid + adsorbed phases).
"""
function _total_column_moles(state, V_fluid, V_solid, comp)
    P = state[:Pressure]
    T = state[:Temperature]
    y = state[:y]
    q = state[:AdsorbedConcentration]
    nc = length(P)
    n = 0.0
    for c in 1:nc
        n += y[comp, c] * P[c] / (GAS_CONSTANT * T[c]) * V_fluid[c]
        n += q[comp, c] * V_solid[c]
    end
    return n
end

"""
Compute thermal energies stored in the column at a given state.
Returns `(E_solid, E_fluid, E_adsorbed)` in Joules.

- `E_solid`:    Σ ρ_s · C_ps · T · V_solid
- `E_fluid`:    Σ C̄_pg · M̄ · P / R · V_fluid  (= C̄_pg · M̄ · C_tot · T · V_fluid)
- `E_adsorbed`: Σ C̄_pa · M̄ · q_tot · T · V_solid
"""
function _column_thermal_energies(state, model::AdsorptionModel, V_fluid, V_solid)
    P = state[:Pressure]
    T = state[:Temperature]
    y = state[:y]
    q = state[:AdsorbedConcentration]
    sys = model.system
    mm = sys.molecular_masses
    C_pg_sys = sys.heat_capacity_gas
    C_pa_sys = sys.heat_capacity_adsorbed
    nc = length(P)
    N = number_of_components(sys)
    ρ_s = first(model.data_domain[:adsorbent_density, Column()])
    C_ps = first(model.data_domain[:adsorbent_heat_capacity, Column()])

    E_solid = 0.0
    E_fluid = 0.0
    E_ads   = 0.0

    for c in 1:nc
        T_c = T[c]
        P_c = P[c]
        C_pg_c = sum(y[i, c] * C_pg_sys[i] for i in 1:N)
        avm_c  = sum(y[i, c] * mm[i]        for i in 1:N)
        C_pa_c = sum(y[i, c] * C_pa_sys[i]  for i in 1:N)
        q_tot_c = sum(q[i, c] for i in 1:N)

        E_solid += ρ_s * C_ps * T_c * V_solid[c]
        E_fluid += C_pg_c * avm_c * P_c / GAS_CONSTANT * V_fluid[c]
        E_ads   += C_pa_c * avm_c * q_tot_c * T_c * V_solid[c]
    end
    return (E_solid, E_fluid, E_ads)
end

# --- Boundary molar flux rates (mol/s) ---
# Returns `(left_net, right_net)` vectors of length N.
# Positive value = net molar flow **into** the column; negative = out of column.

function _boundary_molar_rates(bc::AdsorptionBC, state, model::AdsorptionModel, time)
    N   = number_of_components(model.system)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    nc  = Jutul.number_of_cells(model.domain)
    P   = state[:Pressure]
    T   = state[:Temperature]
    y   = state[:y]

    Af          = _column_face_area_from_model(model)
    trans_left  = _bc_trans_from_model(model, 1)
    trans_right = _bc_trans_from_model(model, nc)

    # Left boundary – fixed-velocity inlet
    q_left    = -bc.v_feed * Af
    P_bc_left = q_left / (trans_left * mob) + P[1]
    c_tot_bc  = P_bc_left / (bc.T_feed * GAS_CONSTANT)

    left_net = zeros(N)
    for i in 1:N
        mass_flux_i = c_tot_bc * q_left * (bc.y_feed[i] - y[i, 1]) +
                      q_left * (bc.y_feed[i] * c_tot_bc)
        left_net[i] = -mass_flux_i   # positive = into column
    end

    # Right boundary – pressure outlet at PH
    q_right     = -trans_right * mob * (bc.PH - P[nc])
    c_tot_right = P[nc] / (T[nc] * GAS_CONSTANT)
    right_net   = [-(q_right * y[i, nc] * c_tot_right) for i in 1:N]

    return (left_net, right_net)
end

function _boundary_molar_rates(bc::PressurisationBC, state, model::AdsorptionModel, time)
    N   = number_of_components(model.system)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    P   = state[:Pressure]
    y   = state[:y]

    trans_left = _bc_trans_from_model(model, 1)
    P_bc_left  = pressure_left(bc, time)
    q_left     = -trans_left * mob * (P[1] - P_bc_left)
    c_tot_bc   = P_bc_left / (bc.T_feed * GAS_CONSTANT)

    left_net = zeros(N)
    for i in 1:N
        mass_flux_i = c_tot_bc * q_left * (bc.y_feed[i] - y[i, 1]) +
                      q_left * (bc.y_feed[i] * c_tot_bc)
        left_net[i] = -mass_flux_i
    end

    return (left_net, zeros(N))
end

function _boundary_molar_rates(bc::BlowdownBC, state, model::AdsorptionModel, time)
    N   = number_of_components(model.system)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    nc  = Jutul.number_of_cells(model.domain)
    P   = state[:Pressure]
    T   = state[:Temperature]
    y   = state[:y]

    trans_right  = _bc_trans_from_model(model, nc)
    P_bc_right   = pressure_right(bc, time)
    q_right      = -trans_right * mob * (P_bc_right - P[nc])
    c_tot_right  = P[nc] / (T[nc] * GAS_CONSTANT)
    right_net    = [-(q_right * y[i, nc] * c_tot_right) for i in 1:N]

    return (zeros(N), right_net)
end

function _boundary_molar_rates(bc::EvacuationBC, state, model::AdsorptionModel, time)
    N   = number_of_components(model.system)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    P   = state[:Pressure]
    T   = state[:Temperature]
    y   = state[:y]

    trans_left = _bc_trans_from_model(model, 1)
    P_bc_left  = pressure_left(bc, time)
    q_left     = -trans_left * mob * (P_bc_left - P[1])
    c_tot_left = P[1] / (T[1] * GAS_CONSTANT)
    # mass_flux = -q * y * cTot; acc -= mass_flux → net inflow = mass_flux
    left_net = [-q_left * y[i, 1] * c_tot_left for i in 1:N]

    return (left_net, zeros(N))
end

# --- Boundary energy flux rates (J/s) ---
# Returns `(Q_left, Q_right)` scalars.
# Positive = energy flowing **into** the column; negative = energy leaving.
#
# Derived from `apply_forces_to_equation!` for `:ColumnConservedEnergy`:
#   `acc[cell] -= bc_src`  →  rate of energy added to cell = bc_src

function _boundary_energy_rates(bc::AdsorptionBC, state, model::AdsorptionModel, time)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    nc  = Jutul.number_of_cells(model.domain)
    ρ_g = first(model.data_domain[:fluid_density, Column()])

    P   = state[:Pressure]
    T   = state[:Temperature]
    y   = state[:y]
    sys = model.system
    N   = number_of_components(sys)

    Af          = _column_face_area_from_model(model)
    trans_left  = _bc_trans_from_model(model, 1)
    trans_right = _bc_trans_from_model(model, nc)

    # Left: fixed-velocity inlet
    q_left    = -bc.v_feed * Af
    P_bc_left = q_left / (trans_left * mob) + P[1]
    C_pg_left = sum(y[i, 1] * sys.heat_capacity_gas[i] for i in 1:N)
    avm_left  = sum(y[i, 1] * sys.molecular_masses[i]  for i in 1:N)
    Q_left    = -((q_left * ρ_g * C_pg_left * (bc.T_feed - T[1])) +
                  (q_left * P_bc_left / GAS_CONSTANT) * C_pg_left * avm_left)

    # Right: pressure outlet
    q_right     = -trans_right * mob * (bc.PH - P[nc])
    C_pg_right  = sum(y[i, nc] * sys.heat_capacity_gas[i] for i in 1:N)
    avm_right   = sum(y[i, nc] * sys.molecular_masses[i]  for i in 1:N)
    Q_right     = -(q_right * P[nc] / GAS_CONSTANT * C_pg_right * avm_right)

    return (Q_left, Q_right)
end

function _boundary_energy_rates(bc::PressurisationBC, state, model::AdsorptionModel, time)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    ρ_g = first(model.data_domain[:fluid_density, Column()])

    P   = state[:Pressure]
    T   = state[:Temperature]
    y   = state[:y]
    sys = model.system
    N   = number_of_components(sys)

    trans_left = _bc_trans_from_model(model, 1)
    P_bc_left  = pressure_left(bc, time)
    q_left     = -trans_left * mob * (P[1] - P_bc_left)
    C_pg_left  = sum(y[i, 1] * sys.heat_capacity_gas[i] for i in 1:N)
    avm_left   = sum(y[i, 1] * sys.molecular_masses[i]  for i in 1:N)
    Q_left     = -((q_left * ρ_g * C_pg_left * (bc.T_feed - T[1])) +
                   (q_left * P_bc_left / GAS_CONSTANT) * C_pg_left * avm_left)

    return (Q_left, 0.0)
end

function _boundary_energy_rates(bc::BlowdownBC, state, model::AdsorptionModel, time)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ
    nc  = Jutul.number_of_cells(model.domain)

    P   = state[:Pressure]
    y   = state[:y]
    sys = model.system
    N   = number_of_components(sys)

    trans_right = _bc_trans_from_model(model, nc)
    P_bc_right  = pressure_right(bc, time)
    q_right     = -trans_right * mob * (P_bc_right - P[nc])
    C_pg_right  = sum(y[i, nc] * sys.heat_capacity_gas[i] for i in 1:N)
    avm_right   = sum(y[i, nc] * sys.molecular_masses[i]  for i in 1:N)
    Q_right     = -(q_right * P[nc] / GAS_CONSTANT * C_pg_right * avm_right)

    return (0.0, Q_right)
end

function _boundary_energy_rates(bc::EvacuationBC, state, model::AdsorptionModel, time)
    μ   = first(model.data_domain[:fluid_viscosity, Column()])
    mob = 1.0 / μ

    P   = state[:Pressure]
    y   = state[:y]
    sys = model.system
    N   = number_of_components(sys)

    trans_left = _bc_trans_from_model(model, 1)
    P_bc_left  = pressure_left(bc, time)
    q_left     = -trans_left * mob * (P_bc_left - P[1])
    C_pg_left  = sum(y[i, 1] * sys.heat_capacity_gas[i] for i in 1:N)
    avm_left   = sum(y[i, 1] * sys.molecular_masses[i]  for i in 1:N)
    Q_left     = -(q_left * P[1] / GAS_CONSTANT * C_pg_left * avm_left)

    return (Q_left, 0.0)
end

# --- Public API ---

"""
    mass_balance(states, timesteps, model, forces;
                 state0 = nothing, component = 1,
                 t_start = 0.0, t_end = sum(timesteps))

Compute the dimensionless mass balance error for a single component over a time period.

The balance is defined as:

    ((moles_in - moles_out) - Δn_col) / |Δn_col|

where `moles_in` and `moles_out` are the integrated molar fluxes into and out of the
column at all boundaries, and `Δn_col = n_end - n_start` is the change in total moles
(fluid + adsorbed) within the column.  A value close to 0 indicates good conservation.

# Arguments
- `states`:     Vector of state snapshots from [`simulate_process`](@ref).
- `timesteps`:  Vector of timestep sizes [s] aligned with `states`.
- `model`:      The [`AdsorptionModel`](@ref).
- `forces`:     Forces applied during simulation – either a single `NamedTuple` or a
                vector of `NamedTuple` (one per step), as returned by [`setup_forces`](@ref).

# Keyword arguments
- `state0`:    Initial state before the first timestep.  Required when `t_start = 0`;
               if `nothing` the first entry of `states` is used as the starting point.
- `component`: Index of the component to check (default `1`, typically CO₂).
- `t_start`:   Start of the balance period [s] (default `0.0`).
- `t_end`:     End   of the balance period [s] (default `sum(timesteps)`).
"""
function mass_balance(states, timesteps, model::AdsorptionModel, forces;
    state0    = nothing,
    component = 1,
    t_start   = 0.0,
    t_end     = sum(timesteps),
)
    # Pre-compute volumes from model domain
    porosity = model.data_domain[:porosity]
    volumes  = model.data_domain[:volumes]
    V_fluid  = volumes .* porosity
    V_solid  = volumes .* (1.0 .- porosity)

    # Cumulative output times aligned with states
    times = cumsum(timesteps)

    # Step indices that fall within [t_start, t_end]
    step_start = searchsortedfirst(times, t_start)
    step_end   = searchsortedlast(times, t_end)
    step_start = clamp(step_start, 1, length(states))
    step_end   = clamp(step_end,   1, length(states))

    # State just before the period starts
    state_before = if step_start == 1
        isnothing(state0) ? states[1] : state0
    else
        states[step_start - 1]
    end

    n_start = _total_column_moles(state_before, V_fluid, V_solid, component)
    n_end   = _total_column_moles(states[step_end], V_fluid, V_solid, component)
    Δn_col  = n_end - n_start

    # Absolute time at the start of the integration window
    t0 = step_start == 1 ? 0.0 : times[step_start - 1]

    # Integrate boundary molar fluxes over the period
    n_in  = 0.0
    n_out = 0.0
    elapsed = 0.0
    for i in step_start:step_end
        dt    = timesteps[i]
        state = states[i]
        force = _get_force_at_step(forces, i)
        bc    = force.bc
        time_i = t0 + elapsed + dt   # absolute time at end of step i

        left_net, right_net = _boundary_molar_rates(bc, state, model, time_i)

        net_flux_i = left_net[component] + right_net[component]
        if net_flux_i > 0
            n_in  += net_flux_i * dt
        else
            n_out += (-net_flux_i) * dt
        end
        elapsed += dt
    end

    denom = abs(Δn_col)
    if denom < eps()
        # Column inventory unchanged; report imbalance relative to total throughput
        denom = max(n_in + n_out, eps())
    end
    return ((n_in - n_out) - Δn_col) / denom
end

"""
    energy_balance(states, timesteps, model, forces;
                   state0 = nothing,
                   t_start = 0.0, t_end = sum(timesteps))

Compute the dimensionless energy balance error over a time period.

The balance is defined as:

    (Q_in + Q_gen - (Q_out + ΔE_solid + ΔE_fluid + ΔE_adsorbed)) /
    |Q_out + ΔE_solid + ΔE_fluid + ΔE_adsorbed|

where:

- `Q_in`          integrated enthalpy entering the column through all boundaries,
                  plus heat transferred from the column wall when T_wall > T_column [J].
- `Q_gen`         heat of adsorption released inside the column [J]:
                  ``Q_\\text{gen} = -\\sum_{c,i} \\Delta H_i \\, \\Delta q_{i,c} \\, V_{\\text{solid},c}``
                  (positive for exothermic adsorption, negative for desorption).
- `Q_out`         integrated enthalpy leaving the column through all boundaries,
                  plus heat lost from the column to the wall when T_column > T_wall [J].
- `ΔE_solid`      change in adsorbent thermal energy:
                  ``\\Delta(\\rho_s C_{ps} T V_{\\text{solid}})``.
- `ΔE_fluid`      change in fluid thermal energy:
                  ``\\Delta(\\bar C_{pg}\\bar M P / R \\, V_{\\text{fluid}})``.
- `ΔE_adsorbed`   change in adsorbed-phase sensible energy:
                  ``\\Delta(\\bar C_{pa}\\bar M q_{\\text{tot}} T V_{\\text{solid}})``.

A value close to 0 indicates good energy conservation.

# Arguments
- `states`:    Vector of state snapshots from [`simulate_process`](@ref).
- `timesteps`: Vector of timestep sizes [s] aligned with `states`.
- `model`:     The [`AdsorptionModel`](@ref).
- `forces`:    Forces applied during simulation – either a single `NamedTuple` or a
               vector of `NamedTuple` (one per step), as returned by [`setup_forces`](@ref).

# Keyword arguments
- `state0`:  Initial state before the first timestep.  Required when `t_start = 0`;
             if `nothing` the first entry of `states` is used as the starting point.
- `t_start`: Start of the balance period [s] (default `0.0`).
- `t_end`:   End   of the balance period [s] (default `sum(timesteps)`).
"""
function energy_balance(states, timesteps, model::AdsorptionModel, forces;
    state0  = nothing,
    t_start = 0.0,
    t_end   = sum(timesteps),
)
    # Pre-compute volumes from model domain
    porosity = model.data_domain[:porosity]
    volumes  = model.data_domain[:volumes]
    V_fluid  = volumes .* porosity
    V_solid  = volumes .* (1.0 .- porosity)

    # Wall geometry and heat transfer coefficients
    r_in   = first(model.data_domain[:r_in, Column()])
    dx     = model.data_domain[:dx]
    h_in   = first(model.data_domain[:inner_htc, Column()])
    WallAreaIn = 2π .* r_in .* dx   # lateral inner surface area per cell [m²]

    sys = model.system
    iso = sys.isotherm
    N   = number_of_components(sys)

    # Cumulative output times aligned with states
    times = cumsum(timesteps)

    # Step indices that fall within [t_start, t_end]
    step_start = searchsortedfirst(times, t_start)
    step_end   = searchsortedlast(times, t_end)
    step_start = clamp(step_start, 1, length(states))
    step_end   = clamp(step_end,   1, length(states))

    # State just before the period starts
    state_before = if step_start == 1
        isnothing(state0) ? states[1] : state0
    else
        states[step_start - 1]
    end

    # Column thermal energies at start and end
    (Es_start, Ef_start, Ea_start) =
        _column_thermal_energies(state_before, model, V_fluid, V_solid)
    (Es_end, Ef_end, Ea_end) =
        _column_thermal_energies(states[step_end], model, V_fluid, V_solid)

    ΔE_solid = Es_end - Es_start
    ΔE_fluid = Ef_end - Ef_start
    ΔE_ads   = Ea_end - Ea_start

    # Heat of adsorption generated over the period
    # Q_gen = -Σ_{c,i} ΔH_i * (q_end - q_start)_i * V_solid (positive = exothermic)
    q_start = state_before[:AdsorbedConcentration]
    q_end   = states[step_end][:AdsorbedConcentration]
    # Use mid-state concentrations for ΔH evaluation (average of start and end)
    P_mid = (state_before[:Pressure] .+ states[step_end][:Pressure]) ./ 2
    T_mid = (state_before[:Temperature] .+ states[step_end][:Temperature]) ./ 2
    y_mid = (state_before[:y] .+ states[step_end][:y]) ./ 2
    nc    = length(P_mid)

    Q_gen = 0.0
    for c in 1:nc
        C_mid = SVector{N, Float64}(y_mid[i, c] * P_mid[c] / (GAS_CONSTANT * T_mid[c])
                                    for i in 1:N)
        ΔH_c  = compute_enthalpy(iso, C_mid, T_mid[c])
        for i in 1:N
            Q_gen += (-ΔH_c[i]) * (q_end[i, c] - q_start[i, c]) * V_solid[c]
        end
    end

    # Absolute time at the start of the integration window
    t0 = step_start == 1 ? 0.0 : times[step_start - 1]

    # Integrate boundary energy fluxes and wall–column heat exchange over the period
    Q_in  = 0.0
    Q_out = 0.0
    elapsed = 0.0
    for i in step_start:step_end
        dt    = timesteps[i]
        state = states[i]
        force = _get_force_at_step(forces, i)
        bc    = force.bc
        time_i = t0 + elapsed + dt

        # Fluid boundary enthalpy fluxes
        Q_L, Q_R = _boundary_energy_rates(bc, state, model, time_i)
        Q_bnd = Q_L + Q_R   # positive = net into column

        if Q_bnd > 0
            Q_in  += Q_bnd * dt
        else
            Q_out += (-Q_bnd) * dt
        end

        # Column–wall heat exchange: source_term = WallAreaIn * h_in * (T - T_w)
        # positive = heat flows from column to wall (energy leaving column)
        T_col  = state[:Temperature]
        T_wall = state[:WallTemperature]
        Q_wall_i = sum(WallAreaIn[c] * h_in * (T_col[c] - T_wall[c]) for c in 1:nc)
        if Q_wall_i > 0
            Q_out += Q_wall_i * dt
        else
            Q_in  += (-Q_wall_i) * dt
        end

        elapsed += dt
    end

    denominator = abs(Q_out + ΔE_solid + ΔE_fluid + ΔE_ads)
    if denominator < eps()
        denominator = max(abs(Q_in + Q_gen), eps())
    end
    return (Q_in + Q_gen - (Q_out + ΔE_solid + ΔE_fluid + ΔE_ads)) / denominator
end
