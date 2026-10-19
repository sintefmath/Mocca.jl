# Moving Mocca 0.1.0 constants structs over to the user interface. Included
# after the interface, whose gas, feed, adsorbent and column it builds.

"""
    convert_constants(constants) → NamedTuple

The [`Gas`](@ref), [`Feed`](@ref), [`Adsorbent`](@ref) and [`Column`](@ref)
described by a `HaghpanahConstants` or `adsorptionConstants` struct, and its
operating conditions: `p_high`, `p_intermediate` and `p_low` [Pa], the
superficial feed velocity `v_feed` and the interstitial `velocity` [m/s], the
ramp rate `λ` [1/s] and the `dispersion_velocity` [m/s] the struct used.
"""
function convert_constants(c::ConstantsStruct)
    names = c isa HaghpanahConstants ? ["CO2", "N2"] : collect(String, c.component_names)
    M = c isa HaghpanahConstants ? [c.molecularMassOfCO2, c.molecularMassOfN2] : collect(c.molecular_masses)
    gas = Gas(names; molar_masses = M, heat_capacities = collect(c.C_pg), viscosity = c.fluid_viscosity, density = c.ρ_g, diffusivity = c.D_m)
    # T0 is the initial temperature of the struct, not the reference
    # temperature of the isotherm, which keeps its default
    iso = dual_site_langmuir((n => (qsb = c.qsbi[i], b0 = c.b0[i], ΔUb = c.ΔUbi[i], qsd = c.qsdi[i], d0 = c.d0[i], ΔUd = c.ΔUdi[i]) for (i, n) in enumerate(names))...)
    adsorbent = Adsorbent(""; isotherm = iso, density = c.ρ_s, heat_capacity = c.C_ps, particle_diameter = c.d_p,
        particle_porosity = c.ϵ_p, tortuosity = c.τ, adsorbed_heat_capacity = collect(c.C_pa))
    column = Column(; length = c.L, radius = c.r_in, outer_radius = c.r_out, porosity = c.Φ, axial_conductivity = c.K_z,
        wall = Wall(density = c.ρ_w, heat_capacity = c.C_pw, conductivity = c.K_w),
        inner_heat_transfer = c.h_in, outer_heat_transfer = c.h_out, ambient_temperature = c.T_a)
    feed = Feed(gas; composition = collect(c.y_feed), temperature = c.T_feed, pressure = c.p_high)
    return (gas = gas, feed = feed, adsorbent = adsorbent, column = column, p_high = c.p_high, p_intermediate = c.p_intermediate,
        p_low = c.p_low, v_feed = c.v_feed, velocity = c.v_feed / c.Φ, λ = c.λ, dispersion_velocity = c.V0_inter)
end
