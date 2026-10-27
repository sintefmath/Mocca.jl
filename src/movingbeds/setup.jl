"""
    tube_bundle(; area_density, sorbent_side_htc, fluid_mass_flow,
        fluid_inlet_temperature, fluid_side_htc = 570.0,
        outer_side = 8e-3, inner_diameter = 6e-3,
        tube_density = 4420.0, tube_heat_capacity = 526.0,
        fluid_density = 1000.0, fluid_heat_capacity = 4200.0)

Properties of the tubes of a [`HeatExchanger`](@ref) section, to pass as
`heat_exchanger` to [`setup_moving_bed_model`](@ref).

`area_density` is the outer tube area per volume of the section [m²/m³]. The
tubes have a square outside of side `outer_side` and a round bore of diameter
`inner_diameter`; from these follow the inner area and the volumes of tube wall
and fluid per volume of the section. The defaults are the 3D-printed TiAl6V4
tubes of the MBTSA demonstrator of Mondino et al. (2022, Table 1), with water
as the fluid and the water-side heat transfer coefficient of their Table 3.

`sorbent_side_htc` and `fluid_side_htc` [W/(m²·K)] are the heat transfer
coefficients on the outside and inside of the tubes; see
[`moving_bed_sorbent_htc`](@ref) for the sorbent side. `fluid_mass_flow`
[kg/s] flows from the bottom of the section to the top, entering at
`fluid_inlet_temperature` [K].
"""
function tube_bundle(;
        area_density,
        sorbent_side_htc,
        fluid_mass_flow,
        fluid_inlet_temperature,
        fluid_side_htc = 570.0,
        outer_side = 8e-3,
        inner_diameter = 6e-3,
        tube_density = 4420.0,
        tube_heat_capacity = 526.0,
        fluid_density = 1000.0,
        fluid_heat_capacity = 4200.0,
    )
    # Per metre of tube: outer perimeter, bore perimeter, wall and bore cross-sections
    outer_perimeter = 4 * outer_side
    bore_area = π * inner_diameter^2 / 4
    tube_length_density = area_density / outer_perimeter
    return (
        outer_area_density = area_density,
        inner_area_density = tube_length_density * π * inner_diameter,
        tube_volume_fraction = tube_length_density * (outer_side^2 - bore_area),
        fluid_volume_fraction = tube_length_density * bore_area,
        tube_density = tube_density,
        tube_heat_capacity = tube_heat_capacity,
        fluid_density = fluid_density,
        fluid_heat_capacity = fluid_heat_capacity,
        sorbent_side_htc = sorbent_side_htc,
        fluid_side_htc = fluid_side_htc,
        fluid_mass_flow = fluid_mass_flow,
        fluid_inlet_temperature = fluid_inlet_temperature,
    )
end

"""
    moving_bed_sorbent_htc(velocity; reference_htc = 117.0, reference_velocity = 7.69e-3, exponent = 0.445)

Heat transfer coefficient [W/(m²·K)] between a moving bed of sorbent and the
tubes of a cross-flow heat exchanger, for a sorbent velocity [m/s]. This is the
correlation `Nu = 0.808·Pe^0.445` of Mondino et al. (2022, eq. 11), written as
a power law through their measurement at the highest sorbent velocity (Table 3,
test 9). It fits all their measurements, 69–117 W/(m²·K) for 2.4–7.7 mm/s,
within 5%. Velocities above about 8 mm/s extrapolate the correlation.
"""
function moving_bed_sorbent_htc(velocity; reference_htc = 117.0, reference_velocity = 7.69e-3, exponent = 0.445)
    return reference_htc * (velocity / reference_velocity)^exponent
end

"""
    setup_moving_bed_model(system::MovingBed; ncells, length, radius, porosity,
        particle_diameter, adsorbent_density, adsorbent_heat_capacity,
        sorbent_mass_flow, solid_fraction = 1 - porosity,
        solid_inlet_loading = zeros(N), solid_inlet_temperature = 298.15,
        fluid_viscosity = 1.8e-5,
        dispersion = 1e-4, thermal_conductivity = 0.1,
        permeability = compute_permeability(porosity, particle_diameter),
        particle_porosity = 0.0, solid_inlet_gas_concentration = zeros(N),
        heat_exchanger = nothing, wall = nothing)

Build the `SimulationModel` of one moving bed section: a vertical column of
`length` [m] and `radius` [m] with `ncells` cells from bottom (1) to top.

- `porosity`: gas volume fraction of the bed. `solid_fraction`: volume fraction
  of sorbent particles; less than `1 − porosity` when the section also holds
  structured packing.
- `particle_diameter` [m] sets the permeability by Kozeny–Carman, as for a
  fixed bed. `adsorbent_density` is the particle density [kg/m³].
- `sorbent_mass_flow` [kg/s] moves down through the section, entering at the
  top with `solid_inlet_loading` [mol/m³ of particles, per component] and
  `solid_inlet_temperature` [K].
- `dispersion` [m²/s] and `thermal_conductivity` [W/(m·K)] are the axial
  dispersion coefficient and effective conductivity.
- `particle_porosity`: pore volume fraction of the particles. The gas in the
  pores moves with the sorbent and counts in the gas held in each cell; it
  enters the top of the section with concentrations
  `solid_inlet_gas_concentration` [mol/m³, per component], unless another
  unit feeds the sorbent. 0 (the default) carries no gas with the sorbent.
- `heat_exchanger`: tube properties from [`tube_bundle`](@ref), required for
  the [`HeatExchanger`](@ref) thermal model.
- `wall`: for the [`WithWall`](@ref) thermal model, a `NamedTuple` with `r_out`,
  `density`, `heat_capacity`, `conductivity`, `inner_htc`, `outer_htc` and
  `ambient_temperature`.

All of these become Jutul parameters, so they can be changed per case and
differentiated.
"""
function setup_moving_bed_model(system::MovingBed;
        ncells,
        length,
        radius,
        porosity,
        particle_diameter,
        adsorbent_density,
        adsorbent_heat_capacity,
        sorbent_mass_flow,
        solid_fraction = 1 - porosity,
        solid_inlet_loading = zeros(number_of_components(system)),
        solid_inlet_temperature = 298.15,
        fluid_viscosity = 1.8e-5,
        dispersion = 1e-4,
        thermal_conductivity = 0.1,
        permeability = compute_permeability(porosity, particle_diameter),
        particle_porosity = 0.0,
        solid_inlet_gas_concentration = zeros(number_of_components(system)),
        heat_exchanger = nothing,
        wall = nothing,
    )
    N = number_of_components(system)
    length_q = Base.length(solid_inlet_loading)
    length_q == N || error("solid_inlet_loading has $length_q values for $N components")
    0 < solid_fraction <= 1 - porosity + 1e-12 || error("solid_fraction ($solid_fraction) must be positive and at most 1 − porosity ($(1 - porosity))")
    if has_heat_exchanger(system)
        isnothing(heat_exchanger) && error("A moving bed with the HeatExchanger thermal model needs heat_exchanger = tube_bundle(...)")
    elseif !isnothing(heat_exchanger)
        error("heat_exchanger is given, but the thermal model is $(thermal_model(system)); use HeatExchanger()")
    end
    if has_wall(system)
        isnothing(wall) && error("A moving bed with the WithWall thermal model needs wall properties")
    end

    mesh = column_mesh(ncells, length, radius)
    nc = Jutul.number_of_cells(mesh)
    domain = Jutul.DataDomain(mesh, porosity = fill(porosity, nc), permeability = fill(permeability, nc))
    domain.entities[Unit()] = 1
    domain[:dx] = map(i -> first(Jutul.cell_dims(mesh, i)), 1:nc)
    domain[:solid_fraction] = fill(solid_fraction, nc)
    # The gas around the particles and in their pores
    domain[:fluid_volume] = domain[:volumes] .* (porosity + solid_fraction * particle_porosity)

    domain[:r_in, Unit()] = radius
    domain[:diffusion_coefficient, Unit()] = dispersion
    domain[:thermal_conductivity, Unit()] = thermal_conductivity
    domain[:adsorbent_density, Unit()] = adsorbent_density
    domain[:adsorbent_heat_capacity, Unit()] = adsorbent_heat_capacity
    domain[:fluid_viscosity, Unit()] = fluid_viscosity
    domain[:sorbent_mass_flow, Unit()] = sorbent_mass_flow
    domain[:solid_inlet_loading, Unit()] = reshape(collect(float.(solid_inlet_loading)), N, 1)
    domain[:solid_inlet_temperature, Unit()] = solid_inlet_temperature
    domain[:particle_porosity, Unit()] = particle_porosity
    domain[:solid_inlet_gas_concentration, Unit()] = reshape(collect(float.(solid_inlet_gas_concentration)), N, 1)

    if !isnothing(heat_exchanger)
        for (k, v) in pairs(heat_exchanger)
            domain[Symbol(:hx_, k), Unit()] = v
        end
    end
    if !isnothing(wall)
        domain[:r_out, Unit()] = wall.r_out
        domain[:wall_density, Unit()] = wall.density
        domain[:wall_heat_capacity, Unit()] = wall.heat_capacity
        domain[:wall_conductivity, Unit()] = wall.conductivity
        domain[:inner_htc, Unit()] = wall.inner_htc
        domain[:outer_htc, Unit()] = wall.outer_htc
        domain[:ambient_temperature, Unit()] = wall.ambient_temperature
    end
    return Jutul.SimulationModel(domain, system, general_ad = true)
end

# Parameters come from the data domain, including the volume of sorbent
# particles, which leaves out any packing
Framework.setup_unit_parameters(model::MovingBedModel; kwarg...) = Jutul.setup_parameters(model; kwarg...)

"""
    setup_unit_state(model::MovingBedModel; pressure, temperature, composition,
        wall_temperature = temperature, tube_temperature = temperature,
        fluid_temperature = temperature)

Uniform initial state of a moving bed, with the sorbent at equilibrium with
the gas. A unit with a [`HeatExchanger`](@ref) also starts its tube wall and
heat transfer fluid at `tube_temperature` and `fluid_temperature`.
"""
function Framework.setup_unit_state(model::MovingBedModel; pressure, temperature, composition,
        wall_temperature = temperature, tube_temperature = temperature, fluid_temperature = temperature)
    kw = _gas_state(model, pressure, temperature, composition, wall_temperature)
    if has_heat_exchanger(model.system)
        kw[:TubeTemperature] = tube_temperature
        kw[:HXFluidTemperature] = fluid_temperature
    end
    return setup_process_state(model; kw...)
end
