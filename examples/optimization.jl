## First we load the necessary modules
import Jutul
import JutulDarcy
import Mocca

function setup_case(prm, step_info = missing)
    # TODO: Handle when only some of the parameters are supplied
    v_feed = prm["v_feed"]
    # p_intermediate = prm["p_intermediate"]
    # p_low = prm["p_low"]
    RealT = typeof(v_feed)
    # p_intermediate#::RealT
    # p_low#::RealT

    # We define parameters, and set up the system and domain as in the [Simulate DCB](simulate_DCB.md) example.
    constants = Mocca.HaghpanahConstants{RealT}(v_feed = v_feed)#, p_intermediate = p_intermediate, p_low = p_low)
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

    # t_press = 2
    # t_ads = 2
    # t_blow = 3
    # t_evac = 4

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
    # TODO: Ensure that these are the correct type from the start?
    total_co2_flux_in = 0.0
    total_co2_flux_out = 0.0

    # TODO: Look at total contributions to total_in and total_out from each process stage
    total_pres = 0.0
    total_ads = 0.0
    total_evac = 0.0

    for (step_info, state, force_outer) in zip(step_infos, states, forces)
        step = step_info[:step]
        dt = step_info[:dt]
        time = step_info[:time]

        if time == 220.0
            @show time
            force = force_outer.bc
            #trans = Mocca.calc_bc_trans(model, state)
            #cell_left = force.cell_left

            # CO2 in for Pressurisation
            if force isa Mocca.PressurisationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_pres -= mass_flux[Mocca.CO2INDEX] * dt
                total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
                # @info "PressurisationBC at $step" mass_flux[Mocca.CO2INDEX]
            end

            # TODO: This is the one causing the gradient mismatch between adjoint and numerical
            # TODO: Create a simpler summed objective function that only looks at Adsorption
            # CO2 in for Adsorption
            if force isa Mocca.AdsorptionBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_ads -= mass_flux[Mocca.CO2INDEX] * dt
                # total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
                # @info "Adsorption at $step" mass_flux[Mocca.CO2INDEX]


                cell_left = 1
                pars = model.system.p
                R = pars.R
                mob = 1.0 / pars.fluid_viscosity
                trans = Mocca.calc_bc_trans(model, state)

                P = state.Pressure[cell_left]
                y = state.y[:, cell_left]
                y_bc = force.y_feed
                T_bc = force.T_feed

                q = Mocca.flux_left(model, state, force)
                P_bc = q / (trans * mob) + P

                c_tot = P_bc / (T_bc * R)
                c = y_bc .* c_tot

                # TODO: q is the culprit! (As a result of force.v_feed in flux_left not being AD/having zero gradient??)
                # Does not seem to be included in sparsity tracing, and adjoint gradient wrt q is 0.0
                # The comes from the fact that v_feed is not an AD variable
                #mass_flux = c_tot .* q .* (y_bc .- y) .+ q .* c
                mass_flux = c_tot[1] * q[1] * (y_bc[1] - y[1]) + q[1] * c[1]
                # @info q, P
                # obj_val = q
                # total_co2_flux_in = mass_flux * dt
                # @show q
                total_co2_flux_in = q
            end

            if force isa Mocca.EvacuationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_evac -= mass_flux[Mocca.CO2INDEX] * dt
                total_co2_flux_out -= mass_flux[Mocca.CO2INDEX] * dt
            end
        end
    end

    # recovery = total_co2_flux_out# /(total_co2_flux_in + total_co2_flux_out)
    recovery = total_co2_flux_in
    @show recovery
    # recovery = total_co2_flux_out/total_co2_flux_in
    return recovery
end
wrapped_global_objective = Jutul.WrappedGlobalObjective(global_objective)

function summed_objective(model, state, dt, step_info, force_outer)

    step = step_info[:step]
    dt = step_info[:dt]
    time = step_info[:time]
    obj_val = 0.0
    # obj_val = 0.0*model.system.p.v_feed

    if time >= 200.0
        # @show model.system.p.v_feed
        force = force_outer.bc

        # CO2 in for Pressurisation
        # if force isa Mocca.PressurisationBC
        #     mass_flux = Mocca.mass_flux_left(state, model, time, force)
        #     obj_val = -mass_flux[Mocca.CO2INDEX] * dt
        #     # @info "PressurisationBC at $step" mass_flux[Mocca.CO2INDEX]
        # end

        # TODO: All of the Adsorption steps are consistently wrong (individually)
        # TODO: Check within mass_flux_left to see if something fishy is going on
        # CO2 in for Adsorption
        if force isa Mocca.AdsorptionBC
            @show time
            #mass_flux = Mocca.mass_flux_left(state, model, time, force)
            #obj_val = -mass_flux[Mocca.CO2INDEX] * dt

            cell_left = 1
            pars = model.system.p
            R = pars.R
            mob = 1.0 / pars.fluid_viscosity
            trans = Mocca.calc_bc_trans(model, state)

            P = state.Pressure[cell_left]
            y = state.y[:, cell_left]
            y_bc = force.y_feed
            T_bc = force.T_feed

            q = Mocca.flux_left(model, state, force)
            P_bc = q / (trans * mob) + P

            c_tot = P_bc / (T_bc * R)
            c = y_bc .* c_tot

            # TODO: q is the culprit! (As a result of force.v_feed in flux_left not being AD/having zero gradient??)
            # Does not seem to be included in sparsity tracing, and adjoint gradient wrt q is 0.0
            # The comes from the fact that v_feed is not an AD variable
            #mass_flux = c_tot .* q .* (y_bc .- y) .+ q .* c
            mass_flux = c_tot[1]*q[1]*(y_bc[1]-y[1]) + q[1]*c[1]
            # @info q, P
            # obj_val = q
            obj_val = mass_flux

            #@info "Adsorption at $step/$time:" mass_flux[Mocca.CO2INDEX]
        end

        # if force isa Mocca.EvacuationBC
        #     mass_flux = Mocca.mass_flux_left(state, model, time, force)
        #     obj_val = -mass_flux[Mocca.CO2INDEX] * dt
        # end
    end
    return obj_val
end
wrapped_sum_objective = Jutul.WrappedSumObjective(summed_objective)

constants_ref = Mocca.HaghpanahConstants{Float64}()
prm_guess = Dict(
    "v_feed" => constants_ref.v_feed,
    #"p_intermediate" => constants_ref.p_intermediate,
    #"p_low" => constants_ref.p_low
)

bar = Jutul.si_unit(:bar)
dprm = Jutul.DictOptimization.DictParameters(prm_guess)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "v_feed"; abs_min = 0.1, abs_max = 2.0)
#Jutul.DictOptimization.free_optimization_parameter!(dprm, "p_intermediate"; abs_min = 0.05bar, abs_max = 0.5bar)
#Jutul.DictOptimization.free_optimization_parameter!(dprm, "p_low"; abs_min = 0.05bar, abs_max = 0.5bar)


## Look at adjoint vs numerical gradient
# New version, using JutulOptimizationProblem
opt_problem = Jutul.DictOptimization.JutulOptimizationProblem(dprm, wrapped_global_objective, setup_case)
obj_and_dobj_adj = Jutul.DictOptimization.evaluate(opt_problem)
dobj_finite_diff = Jutul.DictOptimization.finite_difference_gradient_entry(opt_problem)
println("Numerical: $dobj_finite_diff, adjoint: $(only(obj_and_dobj_adj[2]))")

##
# Old version
"""
grad_adj = Jutul.DictOptimization.parameters_gradient(dprm, wrapped_global_objective, setup_case; raw_output = true)

eps = 1e-3

prm_0 = Dict()
prm_0["v_feed"] = 0.37
case_0 = setup_case(prm_0)
result_0 = Jutul.simulate(case_0; output_substates = true, info_level=0)
packed_steps_0 = Jutul.AdjointPackedResult(result_0, case_0.forces)
obj_val_0 = Jutul.evaluate_objective(wrapped_global_objective, case_0.model, packed_steps_0)

prm_1 = Dict()
prm_1["v_feed"] = prm_0["v_feed"] + eps
case_1 = setup_case(prm_1)
result_1 = Jutul.simulate(case_1; output_substates = true, info_level=0)
packed_steps_1 = Jutul.AdjointPackedResult(result_1, case_1.forces)
obj_val_1 = Jutul.evaluate_objective(wrapped_global_objective, case_1.model, packed_steps_1)

grad_num = (obj_val_1-obj_val_0)/eps"""
nothing

##
prm_opt = Jutul.DictOptimization.optimize(dprm, wrapped_global_objective, setup_case;
    max_it=10,
    maximize=true,
    info_level=-1
)
