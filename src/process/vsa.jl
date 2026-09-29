"""
    four_stage_vsa_flowsheet(constants; ncells = 200, system = FixedBed(constants),
        stage_durations = [15.0, 15.0, 30.0, 40.0], conductance = nothing)

The single-bed, four-stage VSA cycle (pressurisation, adsorption, blowdown,
evacuation) of Haghpanah et al. (2013) as a flowsheet:

    feed ──V_feed──┐                 ┌──V_product── raffinate
                   bottom  Bed  top
    extract ──V_vacuum──┘

Returns `(fs, stages)`. Each stage is a set of device settings; nothing inside
the bed changes between stages:

| Stage          | V_feed                        | V_product                 | V_vacuum                  |
|----------------|-------------------------------|---------------------------|---------------------------|
| pressurisation | valve, feed ramps PL → PH     | closed                    | closed                    |
| adsorption     | flow controller, v_feed·A     | valve, raffinate at PH    | closed                    |
| blowdown       | closed                        | valve, raffinate PH → PI  | closed                    |
| evacuation     | closed                        | closed                    | valve, extract PI → PL    |

Valves use `conductance` [m³/(s·Pa)]; by default the half-cell conductance of
the bed ends (see [`port_conductance`](@ref)), which reproduces the boundary
conditions of earlier Mocca versions.

The raffinate and extract lines only matter if gas flows back into the bed. By
default the raffinate holds the feed without its first (most strongly adsorbed)
component, and the extract holds the feed.
"""
function four_stage_vsa_flowsheet(constants;
        ncells = 200,
        system = FixedBed(constants),
        stage_durations = [15.0, 15.0, 30.0, 40.0],
        conductance = nothing,
        raffinate_composition = light_product_composition(constants.y_feed),
        extract_composition = constants.y_feed,
    )
    bed = setup_process_model(system, constants; ncells = ncells)
    names = component_names(system)
    fs = Flowsheet()
    add_unit!(fs, :Bed, bed)
    for d in (:V_feed, :V_product, :V_vacuum)
        add_unit!(fs, d, setup_device_model(names))
    end
    connect!(fs, (:Bed, :bottom) => (:V_feed, :outlet))
    connect!(fs, (:Bed, :top) => (:V_product, :inlet))
    connect!(fs, (:Bed, :bottom) => (:V_vacuum, :inlet))

    y_feed = constants.y_feed
    T_feed = constants.T_feed
    PH, PI, PL, λ = constants.p_high, constants.p_intermediate, constants.p_low, constants.λ
    feed = PortCondition(pressure = PH, temperature = T_feed, composition = y_feed)
    set_boundary!(fs, :V_feed, :inlet, feed)
    set_boundary!(fs, :V_product, :outlet, PortCondition(pressure = PH, temperature = T_feed, composition = raffinate_composition))
    set_boundary!(fs, :V_vacuum, :outlet, PortCondition(pressure = PL, temperature = T_feed, composition = extract_composition))

    nc = Jutul.number_of_cells(bed.domain)
    C_in = isnothing(conductance) ? port_conductance(bed, 1) : conductance
    C_out = isnothing(conductance) ? port_conductance(bed, nc) : conductance
    A = compute_column_face_area(bed, NamedTuple(Jutul.setup_parameters(bed)))

    t_press, t_ads, t_blow, t_evac = stage_durations
    stages = [
        Stage("pressurisation", t_press;
            V_feed = (law = LinearValve(C_in),
                inlet = PortCondition(pressure = ExponentialRamp(PL, PH, λ), temperature = T_feed, composition = y_feed))),
        Stage("adsorption", t_ads;
            V_feed = (law = VolumetricFlow(constants.v_feed * A),),
            V_product = (law = LinearValve(C_out),)),
        Stage("blowdown", t_blow;
            V_product = (law = LinearValve(C_out),
                outlet = PortCondition(pressure = ExponentialRamp(PH, PI, λ), temperature = T_feed, composition = raffinate_composition))),
        Stage("evacuation", t_evac;
            V_vacuum = (law = LinearValve(C_in),
                outlet = PortCondition(pressure = ExponentialRamp(PI, PL, λ), temperature = T_feed, composition = extract_composition))),
    ]
    return (fs, stages)
end

"""
    light_product_composition(y_feed; trace = 1e-10)

The feed composition with its first component reduced to `trace` and the rest
renormalised: a stand-in for the light product of a separation.
"""
function light_product_composition(y_feed; trace = 1e-10)
    rest = y_feed[2:end]
    rest = rest ./ sum(rest) .* (1 - trace)
    return vcat(trace, rest)
end

"""
    two_bed_vsa_flowsheet(constants; ncells = 100, system = FixedBed(constants),
        feed_time = 15.0, pressurisation_time = 15.0, equalisation_time = 10.0,
        equalisation_conductance = 1e-6, conductance = nothing)

Two beds, A and B, running the four-stage VSA cycle half a cycle apart, with
a pressure equalisation valve joining their tops:

    feed ──V_feed_A──┐              ┌──V_product_A── raffinate
                     bottom  A  top ┤
    extract ─V_vacuum_A┘            V_eq
    feed ──V_feed_B──┐              │
                     bottom  B  top ┤
    extract ─V_vacuum_B┘            └──V_product_B── raffinate

The cycle has six stages. After a bed finishes adsorption it equalises with
the other, evacuated bed, which saves feed gas and compression work:

| Stage | Bed A                  | Bed B                  |
|-------|------------------------|------------------------|
| 1     | pressurisation (feed)  | evacuation             |
| 2     | adsorption             | evacuation             |
| 3     | equalisation (gives)   | equalisation (takes)   |
| 4     | evacuation             | pressurisation (feed)  |
| 5     | evacuation             | adsorption             |
| 6     | equalisation (takes)   | equalisation (gives)   |

Returns `(fs, stages)`.
"""
function two_bed_vsa_flowsheet(constants;
        ncells = 100,
        system = FixedBed(constants),
        feed_time = 15.0,
        pressurisation_time = 15.0,
        equalisation_time = 10.0,
        equalisation_conductance = 1e-6,
        conductance = nothing,
        raffinate_composition = light_product_composition(constants.y_feed),
        extract_composition = constants.y_feed,
    )
    names = component_names(system)
    fs = Flowsheet()
    for b in (:A, :B)
        add_unit!(fs, b, setup_process_model(system, constants; ncells = ncells))
    end
    y_feed, T_feed = constants.y_feed, constants.T_feed
    PH, PI, PL, λ = constants.p_high, constants.p_intermediate, constants.p_low, constants.λ
    # After equalising, a bed sits roughly midway between the high and low
    # pressures; pressurisation starts its feed ramp from there
    P_eq = (PH + PL) / 2
    for b in (:A, :B)
        feed, product, vacuum = Symbol(:V_feed_, b), Symbol(:V_product_, b), Symbol(:V_vacuum_, b)
        for d in (feed, product, vacuum)
            add_unit!(fs, d, setup_device_model(names))
        end
        connect!(fs, (b, :bottom) => (feed, :outlet))
        connect!(fs, (b, :top) => (product, :inlet))
        connect!(fs, (b, :bottom) => (vacuum, :inlet))
        set_boundary!(fs, feed, :inlet, PortCondition(pressure = PH, temperature = T_feed, composition = y_feed))
        set_boundary!(fs, product, :outlet, PortCondition(pressure = PH, temperature = T_feed, composition = raffinate_composition))
        set_boundary!(fs, vacuum, :outlet, PortCondition(pressure = PL, temperature = T_feed, composition = extract_composition))
    end
    add_unit!(fs, :V_eq, setup_device_model(names))
    connect!(fs, (:A, :top) => (:V_eq, :inlet))
    connect!(fs, (:B, :top) => (:V_eq, :outlet))

    bed = fs[:A]
    nc = Jutul.number_of_cells(bed.domain)
    C_in = isnothing(conductance) ? port_conductance(bed, 1) : conductance
    C_out = isnothing(conductance) ? port_conductance(bed, nc) : conductance
    A = compute_column_face_area(bed, NamedTuple(Jutul.setup_parameters(bed)))

    press(b) = Symbol(:V_feed_, b) => (law = LinearValve(C_in),
        inlet = PortCondition(pressure = ExponentialRamp(P_eq, PH, λ), temperature = T_feed, composition = y_feed))
    feed(b) = Symbol(:V_feed_, b) => (law = VolumetricFlow(constants.v_feed * A),)
    product(b) = Symbol(:V_product_, b) => (law = LinearValve(C_out),)
    evacuate(b) = Symbol(:V_vacuum_, b) => (law = LinearValve(C_in),
        outlet = PortCondition(pressure = ExponentialRamp(PI, PL, λ), temperature = T_feed, composition = extract_composition))
    # The equalisation valve's flow direction follows the pressure difference
    equalise = :V_eq => (law = LinearValve(equalisation_conductance),)

    stages = [
        Stage("A pressurisation, B evacuation", pressurisation_time, press(:A), evacuate(:B)),
        Stage("A adsorption, B evacuation", feed_time, feed(:A), product(:A), evacuate(:B)),
        Stage("equalisation A to B", equalisation_time, equalise),
        Stage("B pressurisation, A evacuation", pressurisation_time, press(:B), evacuate(:A)),
        Stage("B adsorption, A evacuation", feed_time, feed(:B), product(:B), evacuate(:A)),
        Stage("equalisation B to A", equalisation_time, equalise),
    ]
    return (fs, stages)
end
