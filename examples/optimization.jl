## First we load the necessary modules
import Jutul
import Mocca

function setup_case(prm, step_info = missing)

    # Set up simulation constants from optimization parameters
    # TODO: Handle only a subset of prm beeing free
    symb_dict = Dict(Symbol(k) => v for (k, v) in prm)
    RealT = valtype(symb_dict)
    constants = Mocca.HaghpanahConstants{RealT}(; symb_dict...)

    # We define parameters, and set up the system and domain as in the [Simulate DCB](simulate_DCB.md) example.
    permeability = Mocca.compute_permeability(constants)
    axial_dispersion = Mocca.calc_dispersion(constants)
    system = Mocca.TwoComponentAdsorptionSystem(; permeability=permeability, dispersion=axial_dispersion, p=constants)

    ncells = 200
    dx = sqrt(pi * constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = Mocca.mocca_domain(mesh, system)

    # # Create the model
    # Now we can assemble the model which contains the domain and the system of equations.
    model = Jutul.SimulationModel(domain, system; general_ad=true)
    push!(model.output_variables, :concentrations)
    push!(model.output_variables, :CellDx)

    # # Set up the initial state
    bar = Jutul.si_unit(:bar)
    P_init = 1 * bar
    T_init = 298.15
    Tw_init = constants.T_a

    yCO2 = fill(1e-10, ncells)
    y_init = hcat(yCO2, 1 .- yCO2)

    state0, parameters = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)

    # # Set up the stage timings

    # Here we have 4 stages and we specify a duration in seconds that we will run each stage.

    t_press = 15
    t_ads = 15
    t_blow = 30
    t_evac = 40

    t_stage = [t_press, t_ads, t_blow, t_evac]

    # We also calculate the time taken to run one cycle and the total time at the end of each stage.

    cycle_time = sum(t_stage)
    step_end = cumsum(t_stage)

    # # Set up boundary conditions

    # We use different boundary conditions for each stage as described above.
    # These boundary conditions are already programmed in Mocca. We just need to define their parameters.

    d_press = Mocca.PressurisationBC(y_feed=constants.y_feed, PH=constants.p_high, PL=constants.p_low,
        λ=constants.λ, T_feed=constants.T_feed, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=0)

    d_ads = Mocca.AdsorptionBC(y_feed=constants.y_feed, PH=constants.p_high, v_feed=constants.v_feed,
        T_feed=constants.T_feed, cell_left=1, cell_right=ncells)

    d_blow = Mocca.BlowdownBC(PH=constants.p_high, PI=constants.p_intermediate,
        λ=constants.λ, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=step_end[2])


    d_evac = Mocca.EvacuationBC(PL=constants.p_low, PI=constants.p_intermediate,
        λ=constants.λ, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=step_end[3])


    # We collect the boundary conditions in the order of their associated stages
    bcs = [d_press, d_ads, d_blow, d_evac]

    # Define the full cyclic simulation by stacking subsequent stages in time
    # for a specified number of cycles
    numcycles = 3

    # TODO: Does this need to have element type inferred?
    timesteps = Float64[]
    sim_forces = []
    maxdt = 1

    for j = 1:numcycles
        for i in eachindex(t_stage)
            numsteps = t_stage[i] / maxdt
            append!(timesteps, repeat([maxdt], Int(floor(numsteps))))
            append!(sim_forces, repeat([Jutul.setup_forces(model, bc=bcs[i])], Int(floor(numsteps))))
        end
    end

    return Jutul.JutulCase(model, timesteps, sim_forces; state0 = state0, parameters = parameters)
end

# TODO: Remove
prm = Dict()
prm["v_feed"] = 0.37
prm["p_intermediate"] = 0.2e5
prm["p_low"] = 0.1e5
case = setup_case(prm)
result = Jutul.simulate(case; output_substates = true, info_level=0)
substates, dt, report_index = Jutul.expand_to_ministeps(result);

## Total objective function recovery
function global_objective(model, state0, states, step_infos, forces, input_data)
    total_co2_flux_in = 0.0
    total_co2_flux_out = 0.0

    for (step_info, state, force_outer) in zip(step_infos, states, forces)
        dt = step_info[:dt]
        time = step_info[:time]

        if time >= 200.0
            force = force_outer.bc

            # CO2 in for Pressurisation
            if force isa Mocca.PressurisationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            end

            # CO2 in for Adsorption
            if force isa Mocca.AdsorptionBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            end

            # CO2 out for Evacuation
            if force isa Mocca.EvacuationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_out -= mass_flux[Mocca.CO2INDEX] * dt
            end
        end
    end

    recovery = total_co2_flux_out/total_co2_flux_in
    return recovery
end
wrapped_global_objective = Jutul.WrappedGlobalObjective(global_objective)


constants_ref = Mocca.HaghpanahConstants{Float64}()
prm_guess = Dict(
    "v_feed" => constants_ref.v_feed,
    "p_intermediate" => constants_ref.p_intermediate,
    "p_low" => constants_ref.p_low
)

bar = Jutul.si_unit(:bar)
dprm = Jutul.DictOptimization.DictParameters(prm_guess)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "v_feed"; abs_min = 0.1, abs_max = 2.0)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "p_intermediate"; abs_min = 0.05bar, abs_max = 0.5bar)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "p_low"; abs_min = 0.05bar, abs_max = 0.5bar)


## Optimize
prm_opt = Jutul.DictOptimization.optimize(dprm, wrapped_global_objective, setup_case;
    max_it=10,
    maximize=true,
    info_level=-1
)
