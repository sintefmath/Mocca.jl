# Process metrics computed from device streams and unit inventories.
#
# Streams are measured where they cross a device, which is what the plant
# supplies or receives: the feed through the feed device, products through the
# product and vacuum devices. The composition of a stream is that of the
# upstream device port.

"""
    cycle_window(timesteps, cycle_time, cycle) → UnitRange

Indices of the steps belonging to cycle number `cycle` (1-based) of length
`cycle_time`, for steps of lengths `timesteps` starting at time 0.
"""
function cycle_window(timesteps, cycle_time, cycle)
    t_end = cumsum(timesteps)
    t0, t1 = (cycle - 1) * cycle_time, cycle * cycle_time
    first_step = findfirst(t -> t > t0 + 1e-9, t_end)
    last_step = findlast(t -> t <= t1 + 1e-9, t_end)
    (isnothing(first_step) || isnothing(last_step)) && error("Cycle $cycle is outside the simulated time")
    return first_step:last_step
end

"""
    stream_totals(states, timesteps, device; window = eachindex(states)) → Vector

Moles of each component that crossed `device` from its inlet to its outlet over
the steps in `window`. Components flowing backwards count as negative.
"""
function stream_totals(states, timesteps, device::Symbol; window = eachindex(states))
    N = size(states[first(window)][device][:InletComposition], 1)
    n = zeros(N)
    for k in window
        s = states[k][device]
        F = s[:MolarFlow][1]
        y = F >= 0 ? s[:InletComposition] : s[:OutletComposition]
        for i in 1:N
            n[i] += F * y[i, 1] * timesteps[k]
        end
    end
    return n
end

"Mole fraction of `component` in the stream with component totals `n`."
purity(n; component = 1) = n[component] / sum(n)

"Fraction of `component` in the feed totals `n_feed` that left in the product totals `n_product`."
recovery(n_product, n_feed; component = 1) = n_product[component] / n_feed[component]

"""
    productivity(n_product, duration, sorbent_mass; component = 1)

Moles of `component` produced per kilogram of sorbent per second.
"""
productivity(n_product, duration, sorbent_mass; component = 1) = n_product[component] / (sorbent_mass * duration)

"""
    sorbent_mass(model, parameters)

Mass of sorbent [kg] in a bed.
"""
sorbent_mass(model::SorbentBedModel, parameters) = sum(parameters[:SolidVolume]) * first(parameters[:AdsorbentDensity])

"""
    vacuum_pump_energy(states, timesteps, device; discharge_pressure = 101325.0,
        efficiency = 0.72, γ = 1.4, window = eachindex(states))

Work [J] to compress the gas leaving through `device` from its outlet pressure
up to `discharge_pressure`, for an ideal gas with heat capacity ratio `γ`
compressed adiabatically with isentropic `efficiency`. The gas enters the pump
at the device's inlet temperature. Only steps with forward flow and an outlet
below the discharge pressure count.

The defaults follow the vacuum pump model of Haghpanah et al. (2013).
"""
function vacuum_pump_energy(states, timesteps, device::Symbol;
        discharge_pressure = 101325.0, efficiency = 0.72, γ = 1.4, window = eachindex(states))
    W = 0.0
    k_exp = (γ - 1) / γ
    for k in window
        s = states[k][device]
        F = s[:MolarFlow][1]
        P_suction = s[:OutletPressure][1]
        if F > 0 && P_suction < discharge_pressure
            T = s[:InletTemperature][1]
            W += F * GAS_CONSTANT * T / k_exp * ((discharge_pressure / P_suction)^k_exp - 1) / efficiency * timesteps[k]
        end
    end
    return W
end

"""
    specific_energy_kwh_per_tonne(energy, n_product, molar_mass; component = 1)

Energy per tonne of `component` produced, in kWh/t, given the energy [J], the
product component totals [mol] and the component's molar mass [kg/mol].
"""
function specific_energy_kwh_per_tonne(energy, n_product, molar_mass; component = 1)
    tonnes = n_product[component] * molar_mass / 1000
    return energy / 3.6e6 / tonnes
end

"""
    component_inventory(model, state, parameters) → Vector

Moles of each component held in a sorbent bed, in the gas and on the sorbent.
"""
function component_inventory(model::SorbentBedModel, state, parameters)
    c = state[:y] .* (state[:Pressure] ./ (GAS_CONSTANT .* state[:Temperature]))'
    gas = c * parameters[:FluidVolume]
    adsorbed = state[:AdsorbedConcentration] * parameters[:SolidVolume]
    return vec(gas .+ adsorbed)
end
