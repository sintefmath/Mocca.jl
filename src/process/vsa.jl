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
