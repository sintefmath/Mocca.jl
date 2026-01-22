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

    t_stage, num_cycles = cycle_definition()
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
wrapped_global_objective = Jutul.WrappedGlobalObjective(objective_func, depends_on_parameters = false);

# We use the original parameter values as a starting point for the optimization
constants_ref, info_ref = Mocca.parse_input(Mocca.haghpanah_cyclic_input(); typeT=Float64)

basecase, = Mocca.setup_mocca_case(constants_ref, info_ref)

prm_guess = Dict(
    "v_feed" => constants_ref.v_feed,
    "p_intermediate" => constants_ref.p_intermediate,
    "p_low" => constants_ref.p_low
)

# We create a setup function for making simulation cases.
# This is needed by the optimizer so that it knows how to set up a new simulation
# from the current iteration of the optimization parameters.
function setup_case(prm, step_info = missing; num_cycles = cycle_definition()[2], state0 = basecase.state0)
    # @info "Solving with $num_cycles"
    param_dict_symb = Dict(Symbol(k) => v for (k, v) in prm)
    RealT = valtype(param_dict_symb)
    constants, info = Mocca.parse_input(Mocca.haghpanah_cyclic_input(); typeT=RealT)
    info.num_cycles = num_cycles;
    for (k, v) in param_dict_symb
        print(k)
        print(v)
        setproperty!(constants, Symbol(k), v)
    end
    if false
        case,  = Mocca.setup_mocca_case(constants, info)
    else
        (; model, parameters) = basecase
        if false
            tmp, = Mocca.setup_mocca_case(constants, info)
            dt = tmp.dt
            forces = tmp.forces
        else
            forces, dt = Mocca.setup_forces(model,
                info.stage_durations,
                info.stage_types;
                num_cycles = info.num_cycles,
                # max_dt = info.maxdt,
                max_dt = 1.0,
                constants = constants # !! Use updated constants
            );
        end
        case = Mocca.MoccaCase(model, dt, forces, state0 = state0, parameters = parameters)
    end
    return case
end;

# c = setup_case(prm_guess);
##
using Statistics
# Specify which parameters we wish to optimize and set limits for their final values. Relative change limits can also be specified.
bar = Jutul.si_unit(:bar)

sim, cfg = Mocca.setup_process_simulator(basecase.model, basecase.state0, basecase.parameters,
    info_level = -1, end_report = false);

num_cycles_outer = 10
num_cycles_optimizer = 3
current_parameters = prm_guess
maxit = 10
num_outer_it = 3

dict_opt = missing
for outer_it in 1:num_outer_it
    case_full = setup_case(current_parameters; num_cycles = num_cycles_outer)
    states, = Mocca.simulate_process(case_full, simulator = sim, config = cfg)

    println("Results after $num_cycles_outer cycles:")
    for (k, v) in states[end]
        @info "$k => Mean: $(Statistics.mean(v)) Min: $(minimum(v)) Max: $(maximum(v))"
    end
    setup_case_inner = (arg...) -> setup_case(arg...; state0 = states[end], num_cycles = num_cycles_optimizer)

    dict_opt = DictParameters(current_parameters)
    free_optimization_parameter!(dict_opt, "v_feed"; abs_min = 0.1, abs_max = 2.0)
    free_optimization_parameter!(dict_opt, "p_intermediate"; abs_min = 0.05bar, abs_max = 0.5bar)
    free_optimization_parameter!(dict_opt, "p_low"; abs_min = 0.05bar, abs_max = 0.5bar)

    println("Starting optimization with $maxit iterations:")
    current_parameters = optimize(dict_opt, wrapped_global_objective, setup_case;
        max_it=maxit,
        maximize=true,
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
end

# # Run the optimization

# We call the optimizer provided by Jutul.
# Note that we are maximizing the objective function.
# import Jutul: Simulator, simulator_config
# sim = Simulator(basecase)
# lsolve = LUSolver()
# lsolve = Jutul.GenericKrylov(:bicgstab, preconditioner = Jutul.ILUZeroPreconditioner())
# cfg = simulator_config(sim, linear_solver = lsolve, end_report = false)
timestep_selector_cfg = (
    y=0.05,
    Temperature=10.0,
    Pressure=10.0
)
sim, cfg = Mocca.setup_process_simulator(basecase.model, basecase.state0, basecase.parameters,
    timestep_selector_cfg = timestep_selector_cfg,
    output_substates = true,
    info_level = 2
);
##
prm_opt = optimize(dprm, wrapped_global_objective, setup_case;
    max_it=10,
    maximize=true,
    # info_level=-1,
        backend_arg = (
        use_sparsity = true,
        di_sparse = true,
        single_step_sparsity = :unique_forces,
        do_prep = true,
        deps = :case,
        deps_ad = :di
    ),
    # simulator = sim,
    # config = cfg
)

# We can plot the optimization history to see how the objective function has changed throughout the optimization
Mocca.plot_optimization_history(dprm; yscale = identity, ylabel = "Recovery")

# Finally, we look at the optimized parameters.
# We see that the optimized intermediate and low pressure values have reached their prescribed limits,
# meaning that we could have increased the objective function further if we were allowed to change the limits.
dprm


# Optimization: Objective #12: 3.95121e-01 (f/f0=1.870e+00), gradient 2-norm: 1.21684e-04
#   10 | 1.7712e+01 | 1.5777e-01 | 1
# Optimization: Finished in 143.1014966 seconds.     
# DictParameters with 3 parameters (3 active), and 0 multipliers:
# Active optimization parameters
# ┌────────────────┬───────────────┬───────┬────────┬─────────┬─────────────────┬────────┐
# │           Name │ Initial value │ Count │    Min │     Max │ Optimized value │ Change │
# ├────────────────┼───────────────┼───────┼────────┼─────────┼─────────────────┼────────┤
# │         v_feed │ 0.37          │     1 │    0.1 │     2.0 │ 0.684           │  85.0% │
# │ p_intermediate │ 20000.0       │     1 │ 5000.0 │ 50000.0 │ 50000.0         │ 150.0% │
# │          p_low │ 10000.0       │     1 │ 5000.0 │ 50000.0 │ 5000.0          │ -50.0% │
# └────────────────┴───────────────┴───────┴────────┴─────────┴─────────────────┴────────┘
# 143 -> 48
##
