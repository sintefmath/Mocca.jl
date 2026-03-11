# # Optimization of simulation models

# We demonstrate optimization of a cyclic Vacuum Swing Adsorption (VSA) simulation model in Mocca.
# The objective is to maximize the CO2 recovery of the process, i.e., the fraction of the input CO2 we are able to capture.
# We do this by tuning:
# * Feed velocity
# * Low pressure
# * Intermediate pressure
#
# We leverage the powerful and flexible optimization functionality of Jutul.
# For more details about the VSA modelling, see the [Simulate cyclic](simulate_cyclic.md) example.

# # Setting up the optimization problem

# Start by importing the necessary modules
import Mocca
import Jutul
import Jutul.DictOptimization: optimize, DictParameters, free_optimization_parameter!


# Create a helper function for getting timing for the stages and the number of cycles
function cycle_definition()
    t_press = 15.0
    t_ads = 15.0
    t_blow = 30.0
    t_evac = 40.0
    t_stage = [t_press, t_ads, t_blow, t_evac]
    num_cycles = 3
    # num_cycles = 500
    return (t_stage, num_cycles)
end;

# Define the objective function. We need access to all timesteps at the same time to calculate the recovery.
# Jutul allows us to do this using a global objective function.
function objective_func(model, state0, states, step_infos, forces, input_data)
    total_co2_flux_in = 0.0
    total_co2_flux_out = 0.0

    t_stage, = cycle_definition()
    num_cycles = model.extra[:num_cycles]

    start_time_last_cycle = sum(t_stage)*(num_cycles-1)

    for (step_info, state, force_outer) in zip(step_infos, states, forces)
        dt = step_info[:dt]
        time = step_info[:time]

        force = force_outer.bc

        co2_flux_in = 0.0
        co2_flux_out = 0.0

        if force isa Mocca.PressurisationBC
            mass_flux = Mocca.mass_flux_left(state, model, time, force)
            # total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            co2_flux_in = mass_flux[Mocca.CO2INDEX] * dt
        end

        if force isa Mocca.AdsorptionBC
            mass_flux = Mocca.mass_flux_left(state, model, time, force)
            # total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            co2_flux_in = mass_flux[Mocca.CO2INDEX] * dt
        end

        if force isa Mocca.EvacuationBC
            mass_flux = Mocca.mass_flux_left(state, model, time, force)
            # total_co2_flux_out -= mass_flux[Mocca.CO2INDEX] * dt
            co2_flux_out = mass_flux[Mocca.CO2INDEX] * dt
        end

        if time >= start_time_last_cycle # We only use the last cycle for calculating the objective, once the system has more or less stabilized
            w = 1.0
        else
            w = 1e-12
        end
        total_co2_flux_out -= w*co2_flux_out
        total_co2_flux_in -= w*co2_flux_in
    end

    recovery = total_co2_flux_out/total_co2_flux_in
    return recovery
end


# We create a setup function for making simulation cases.
# This is needed by the optimizer so that it knows how to set up a new simulation
# from the current iteration of the optimization parameters.
function setup_case(prm, step_info = missing; num_cycles = cycle_definition()[2], state0 = basecase.state0, absLims = missing)
    # @info "Solving with $num_cycles"
    param_dict_symb = Dict(Symbol(k) => v for (k, v) in prm)
    paramLims_dict_symb = Dict(Symbol(k) => v for (k, v) in absLims)

    RealT = valtype(param_dict_symb)
    constants, info = Mocca.parse_input(Mocca.haghpanah_cyclic_input(); typeT=RealT)
    for (k, v) in param_dict_symb
        print(k)
        print(v)
    
        if !ismissing(absLims)
            if v >= paramLims_dict_symb[k][2]
                v = paramLims_dict_symb[k][2]
            elseif v <= paramLims_dict_symb[k][1]
                v = paramLims_dict_symb[k][1]
            end
        end
        setproperty!(constants, Symbol(k), v)
    end


    (; model, parameters) = basecase

        forces, dt = Mocca.setup_forces(model,
            info.stage_durations,
            info.stage_types;
            num_cycles = num_cycles,
            # max_dt = info.maxdt,
            max_dt = 1.0,
            constants = constants # !! Use updated constants
        );

    case = Mocca.MoccaCase(model, dt, forces, state0 = state0, parameters = parameters)

    case.model.extra[:num_cycles] = num_cycles
    return case
end;


function find_steady_state(case, sim, cfg)
    cfg[:info_level] = 0
    states, = Mocca.simulate_process(case, simulator = sim, config = cfg)
    restart = Mocca.setup_process_state(case.model; state0 = states[end])
    return restart
end

# c = setup_case(prm_guess);
##
using Statistics

wrapped_global_objective = Jutul.WrappedGlobalObjective(objective_func, depends_on_parameters = false);

# We use the original parameter values as a starting point for the optimization
constants_ref, info_ref = Mocca.parse_input(Mocca.haghpanah_cyclic_input(); typeT=Float64)

basecase, = Mocca.setup_mocca_case(constants_ref, info_ref)


# Specify which parameters we wish to optimize and set limits for their final values. Relative change limits can also be specified.
bar = Jutul.si_unit(:bar)

prm_guess = Dict(
    "v_feed" => constants_ref.v_feed,
    "p_intermediate" => constants_ref.p_intermediate,
    "p_low" => constants_ref.p_low
)

prm_absLims = Dict(
    "v_feed" => [0.1, 2.0],
    "p_intermediate" => [0.05bar, 0.5bar],
    "p_low" => [0.05bar, 0.5bar]
)



sim, cfg = Mocca.setup_process_simulator(basecase.model, basecase.state0, basecase.parameters,
    info_level = -1, end_report = false);

num_cycles_outer = 2
num_cycles_optimizer = 4
current_parameters = prm_guess
maxit = 10
num_outer_it = 3

dict_opt = missing
results = []
push!(results, current_parameters)
histories = []
initial_states = []
A = [0.0 1.0 -1.0]
b = [0.01bar]
lin_ineq = (A = A, b = b)  # interpreted as Au \leq b

for outer_it in 1:num_outer_it

    restart = find_steady_state(setup_case(current_parameters; num_cycles = num_cycles_outer, absLims = prm_absLims), sim, cfg)
    setup_case_inner = (arg...) -> setup_case(arg...; state0 = restart, num_cycles = num_cycles_optimizer, absLims = prm_absLims)

    global dict_opt = DictParameters(current_parameters, verbose = false)
    free_optimization_parameter!(dict_opt, "v_feed"; abs_min = prm_absLims["v_feed"][1], abs_max = prm_absLims["v_feed"][2])
    free_optimization_parameter!(dict_opt, "p_intermediate"; abs_min = prm_absLims["p_intermediate"][1], abs_max = prm_absLims["p_intermediate"][2])
    free_optimization_parameter!(dict_opt, "p_low"; abs_min = prm_absLims["p_low"][1], abs_max = prm_absLims["p_low"][2])

    println("Starting optimization with $maxit iterations:")
    global current_parameters = optimize(dict_opt, wrapped_global_objective, setup_case_inner;
        max_it=maxit,
        maximize=true,
        lin_ineq = lin_ineq,
        info_level=-1,
            backend_arg = (
            use_sparsity = true,
            di_sparse = true,
            single_step_sparsity = :unique_forces,
            do_prep = true,
            deps = :case,
            deps_ad = :di
        ),
        simulator = sim,
        config = cfg
    )
    push!(results, current_parameters)
    push!(histories, dict_opt.history)
    push!(initial_states, restart)
end

##
using CairoMakie
offset = 0
fig = Figure(size = (800, 800))
ax = Axis(fig[1, 1];
    xlabel = "Iteration",
    ylabel = "Objective",
    title = "Optimization history",
    yscale = log10
)
for (i, history) in enumerate(histories)
    iterations = 1:length(history.objectives)
    objectives = history.objectives
    scatter!(ax, iterations .+ offset, objectives; label = "Outer it $i")
    global offset += length(iterations)
end
axislegend(position = :rb)
ax = Axis(fig[2, 1];
    xlabel = "Parameter",
    ylabel = "Value",
    title = "Optimized parameters"
)
relstate0(k, offset = 0) = map(x -> mean(x[k] .- offset), initial_states)/mean(initial_states[1][k] .- offset)

scatterlines!(ax, relstate0(:Pressure), label = "Pressure")
scatterlines!(ax, relstate0(:Temperature, 273.15), label = "Temperature")
scatterlines!(ax, relstate0(:y), label = "y")
scatterlines!(ax, relstate0(:AdsorbedConcentration), label = "AdsorbedConcentration")

axislegend(position = :lb)
fig
