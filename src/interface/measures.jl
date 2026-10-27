# Performance measures of a process as functions of integrals over its cycle,
# or rates at steady state, and of quantities that do not depend on the state.
# performance evaluates them, and performance_gradient differentiates them by
# the chain rule through the integrals (see measure_partials).

# Units with a heat exchanger whose fluid comes from outside the process, those
# whose fluid comes from another unit only pass heat around. The fluid of a
# hot utility enters above the ambient temperature, that of a cold one at or
# below it.
function _utility_units(fs, parameters, o)
    external = [u for u in holdup_units(fs) if has_heat_exchanger(fs[u].system) && isnothing(Operation._hx_upstream(fs, u))]
    hot = [u for u in external if first(parameters[u][:HXFluidInletTemperature]) > o.ambient_temperature]
    return (hot = hot, cold = setdiff(external, hot))
end

_has_heat_exchangers(fs) = any(u -> has_heat_exchanger(fs[u].system), holdup_units(fs))

"""
    performance_integrals(p) → NamedTuple

The integrals that the performance of process `p` is measured from, as Jutul
sum objectives (see [`StreamObjective`](@ref)), over a cycle, or over one
second at steady state:

- `key_extract`: moles of the key component through the `:extract` devices,
  `total_extract`, of every component through them, and `key_feed`, of the key
  component through the `:feed` devices;
- `energy` [J]: vacuum pumps on the outlet devices whose gas leaves below
  `ambient_pressure`, compressors on the inlet devices whose gas enters above
  it (see [`VacuumPumpObjective`](@ref) and [`CompressorObjective`](@ref)),
  and the fans and blowers of the `:link` devices that move gas up in pressure
  (see [`BlowerObjective`](@ref)); `blower_energy`, those of the links alone. A
  link that lets the gas down, such as the valve of an equalisation, takes no
  work;
- `heating` [J]: heat given by heat transfer fluids that come from outside the
  process hotter than `ambient_temperature`, and `cooling` [J], heat taken by
  those that come in at or below it (see [`HeatDutyObjective`](@ref)); units
  whose fluid comes from another unit only pass heat around;
- `equivalent_work` [J]: the work equivalent to the heating, each bit of heat
  weighted by the Carnot factor `1 − T₀/T` at the temperature `T` of the fluid
  giving it up, against `T₀ = ambient_temperature` (see
  [`EquivalentWorkObjective`](@ref)).

followed by those added with [`add_integral!`](@ref).
"""
function performance_integrals(p::Process)
    fs, o = p.flowsheet, p.options
    k = o.key_component
    N = number_of_components(fs.gas)
    extract = devices_with_role(fs, :extract)
    feed = devices_with_role(fs, :feed)
    energy = Any[]
    # A heat capacity ratio from the feed changes with the feed temperature only
    # if the heat capacities do; that dependence is not differentiated
    γ = pump_heat_capacity_ratio(p)
    for role in OUTLET_ROLES, d in devices_with_role(fs, role)
        push!(energy, VacuumPumpObjective(d; discharge_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = γ))
    end
    for role in INLET_ROLES, d in devices_with_role(fs, role)
        push!(energy, CompressorObjective(d; suction_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = γ))
    end
    blowers = Tuple(BlowerObjective(d; efficiency = o.pump_efficiency, γ = γ) for d in devices_with_role(fs, :link))
    hx = _utility_units(fs, p.parameters, o)
    builtin = (key_extract = SumObjective(Tuple(StreamObjective(d; component = k) for d in extract)),
        total_extract = SumObjective(Tuple(StreamObjective(d; component = i) for d in extract for i in 1:N)),
        key_feed = SumObjective(Tuple(StreamObjective(d; component = k) for d in feed)),
        energy = SumObjective((energy..., blowers...)),
        blower_energy = SumObjective(blowers),
        heating = SumObjective(Tuple(HeatDutyObjective(fs, u) for u in hx.hot)),
        cooling = SumObjective(Tuple(HeatDutyObjective(fs, u; scale = -1.0) for u in hx.cold)),
        equivalent_work = SumObjective(Tuple(EquivalentWorkObjective(fs, u; ambient_temperature = o.ambient_temperature) for u in hx.hot)))
    return merge(builtin, (; p.integrals...))
end

"""
    performance_quantities(p) → NamedTuple

The quantities that the performance of process `p` depends on besides its
integrals (see [`performance_integrals`](@ref)), which do not depend on the
state: the `duration` [s] the integrals are taken over, the cycle time, or 1 s
for a [`SteadyProcess`](@ref), whose integrals are rates; the
`adsorbent_volume` [m³] and `adsorbent_mass` [kg] of its beds; and the
`sorbent_inventory` [kg] of its units whose sorbent moves.
"""
function performance_quantities(p::CyclicProcess)
    return (duration = cycle_time(p), adsorbent_volume = adsorbent_volume(p), adsorbent_mass = adsorbent_mass(p),
        sorbent_inventory = sorbent_inventory(p.flowsheet, p.parameters))
end

function performance_quantities(p::SteadyProcess)
    fs, prm = p.flowsheet, p.parameters
    return (duration = 1.0, adsorbent_volume = sum(b -> adsorbent_volume(fs[b], prm[b]), beds(fs); init = 0.0),
        adsorbent_mass = sum(b -> sorbent_mass(fs[b], prm[b]), beds(fs); init = 0.0), sorbent_inventory = sorbent_inventory(fs, prm))
end

"""
    performance_measures(p) → NamedTuple

The measures of the performance of process `p`, each a function of `s`, the
values of its integrals (see [`performance_integrals`](@ref)) and quantities
(see [`performance_quantities`](@ref)) in one `NamedTuple`:

- `purity` = `s.key_extract / s.total_extract`, mole fraction of the key
  component in the extract;
- `recovery` = `s.key_extract / s.key_feed`;
- `specific_energy` [kWh per tonne of key component in the extract], from
  `s.energy`;
- `productivity` [mol of key component per m³ of adsorbent per s];
- with heat exchangers, `specific_heat_duty` [J per kg of key component], from
  `s.heating`, and `specific_equivalent_work` [kWh per tonne], from
  `s.energy + s.equivalent_work`;

followed by those added with [`add_measure!`](@ref).
"""
function performance_measures(p::Process)
    fs = p.flowsheet
    M = fs.gas.molar_masses[p.options.key_component]
    builtin = (purity = s -> s.key_extract / s.total_extract,
        recovery = s -> s.key_extract / s.key_feed,
        specific_energy = s -> s.energy / 3.6e6 / (s.key_extract * M / 1000),
        productivity = s -> s.key_extract / (s.adsorbent_volume * s.duration))
    if _has_heat_exchangers(fs)
        builtin = merge(builtin, (specific_heat_duty = s -> s.heating / (s.key_extract * M),
            specific_equivalent_work = s -> (s.energy + s.equivalent_work) / 3.6e6 / (s.key_extract * M / 1000)))
    end
    return merge(builtin, (; p.measures...))
end

"""
    performance_sums(p, states, timesteps; window = eachindex(states)) → NamedTuple

The values of the integrals of process `p` (see
[`performance_integrals`](@ref)) over the steps in `window`, followed by its
quantities (see [`performance_quantities`](@ref)), with the `duration` that
of the steps: the argument of each of its measures.
"""
function performance_sums(p::Process, states, timesteps; window = eachindex(states))
    I = integral_values(performance_integrals(p), p.model, states, timesteps, p.parameters; window = window)
    return merge(I, performance_quantities(p), (duration = sum(timesteps[window]),))
end

_measure_values(p, s) = map(f -> f(s), performance_measures(p))

"""
    performance(p::Process; kwarg...) → NamedTuple

Performance of process `p` at the state it settles to, its cyclic steady
state or its steady state; keywords go to [`steady_state`](@ref).
"""
performance(p::Process; kwarg...) = performance(steady_state(p; kwarg...))

"""
    add_integral!(p, name, G)

Add an integral `name` to the performance of process `p`: a Jutul sum
objective `G(model, state, dt, step_info, forces)`, summed over the steps of a
cycle (or taken over one second at steady state), such as
`StreamObjective(:V_product; component = 2)` or `HeatDutyObjective(:Desorption)`.
The state holds each unit's parameters; the objective gets no forces when the
performance is evaluated. Measures added with [`add_measure!`](@ref) can then
use it as `s.name`, and get its gradient from the adjoints.
"""
function add_integral!(p::Process, name::Symbol, G)
    _check_new_name(p, name, "integral")
    push!(p.integrals, name => G)
    return p
end

"""
    add_measure!(p, name, f)

Add a performance measure `name` to process `p`: `f(s)`, a function of the
values `s` of its integrals and quantities (see
[`performance_integrals`](@ref), [`performance_quantities`](@ref) and
[`add_integral!`](@ref)), such as

```julia
add_measure!(p, :nitrogen_in_extract, s -> s.total_extract - s.key_extract)
```

[`performance`](@ref) reports it, and [`performance_gradient`](@ref) and
[`optimise!`](@ref) differentiate it with no derivative code: the chain rule
through `f` is taken by automatic differentiation, so `f` must be a generic
Julia function of numbers.
"""
function add_measure!(p::Process, name::Symbol, f)
    _check_new_name(p, name, "measure")
    push!(p.measures, name => f)
    return p
end

function _check_new_name(p, name, what)
    taken = (keys(performance_integrals(p))..., keys(performance_quantities(p))..., keys(performance_measures(p))...)
    name in taken && input_error("$what $name", "", "the process already has an integral, quantity or measure of that name")
end
