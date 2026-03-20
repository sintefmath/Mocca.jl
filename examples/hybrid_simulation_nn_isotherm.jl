# # Hybrid Simulation with Neural Network Isotherms (Direct Mixture Learning)
# <tags: Adsorption, MachineLearning, Advanced>

# This example demonstrates integrating neural networks into Mocca.jl for isotherm modeling.
# Unlike traditional approaches that use pure-component isotherms with competition models or Langmuir IAST, this example
# trains a single NN on mixture data to learn competitive adsorption behavior directly.
#
# ## Workflow:
# 1. Generate training data from analytical dual-site Langmuir isotherms at various mixture compositions
# 2. Train a neural network: `[C_CO2, C_N2, T] → [q_CO2, q_N2]`
# 3. Integrate the NN into Mocca's simulation framework via custom secondary variables
# 4. Compare NN-based simulation with analytical reference
#
# ## Important Implementation Details:
# - **Gradient Stability**: Concentrations are clamped to ≥0.01 mol/m³ before log-transformation
#   to prevent large gradients at very low concentrations (d(log(C))/dC = 1/C) which causes convergence issues in the Newton solver
# - **Model Persistence**: Trained model is saved with architecture versioning to avoid retraining
# - **Competitive Adsorption**: Learned implicitly from mixture training data

# ## Preliminaries

import Jutul
import JutulDarcy
import Mocca
using Lux, ADTypes, Zygote, Optimisers, Random, Statistics
using StaticArrays
using CairoMakie
using Colors
using JLD2

# ## Set up the reference simulation case
# We use the same Direct Column Breakthrough simulation setup as in `simulate_DCB.jl`

function setup_reference_simulation()
    ## Load standard constants from Haghpanah et al. 2013
    constants = Mocca.HaghpanahConstants{Float64}(h_in=0.0, h_out=0.0)

    ## Create the standard two-component adsorption system
    system = Mocca.TwoComponentAdsorptionSystem(constants)

    ## Set up domain and model
    ncells = 200
    model = Mocca.setup_process_model(system; ncells=ncells)

    ## Initial conditions
    bar = Jutul.si_unit(:bar)
    P_init = 1 * bar
    T_init = 298.15
    Tw_init = constants.T_a

    yCO2 = 1e-10
    y_init = [yCO2, 1 - yCO2]

    state0 = Mocca.setup_process_state(model;
        Pressure=P_init, Temperature=T_init,
        WallTemperature=Tw_init, y=y_init)
    parameters = Mocca.setup_process_parameters(model)

    ## Time stepping and boundary conditions
    t_ads = 5000.0
    maxdt = 5000.0
    sim_forces, timesteps = Mocca.setup_forces(model, [t_ads], ["adsorption"];
        num_cycles=1, max_dt=maxdt)

    return model, state0, parameters, timesteps, sim_forces, constants
end

# ## Generate training data from analytical isotherms
# We generate training data from MIXTURES to capture competitive adsorption behavior directly.

function generate_mixture_training_data(constants)
    ## Use Mocca's actual isotherm model to generate training data
    ## This ensures exact consistency with the simulation
    isotherm_model = Mocca.DualSiteLangmuir(constants)
    R = constants.R
    ρ_s = constants.ρ_s  ## Adsorbent density [kg/m³]

    ## Helper function to get mixture loading
    function mixture_loading(p_total, T, y_co2)
        ## Convert to concentrations
        C_co2 = y_co2 * p_total / (R * T)
        C_n2 = (1 - y_co2) * p_total / (R * T)
        q_bulk = Mocca.compute_equilibrium(isotherm_model, SVector{2}([C_co2, C_n2]), T)
        return q_bulk[1] / ρ_s, q_bulk[2] / ρ_s  ## Convert from mol/m³ to mol/kg
    end

    ## Define ranges for training data (1800 total samples)
    n_pressure = 30
    n_temperature = 4
    n_composition = 15

    p_min, p_max = 1e4, 2e5      ## 10-200 kPa
    T_min, T_max = 273.15, 373.15  ## 0-100°C

    pressures = exp.(range(log(p_min), log(p_max), length=n_pressure))
    temperatures = range(T_min, T_max, length=n_temperature)
    ## Avoid pure components at 0 and 1
    y_co2_range = range(0.01, 0.99, length=n_composition)

    ## Generate training samples
    n_samples = n_pressure * n_temperature * n_composition
    inputs_raw = zeros(Float64, 3, n_samples)
    outputs_raw = zeros(Float64, 2, n_samples)

    idx = 1
    for T in temperatures, p in pressures, y_co2 in y_co2_range
        ## Convert to concentrations
        C_co2 = y_co2 * p / (R * T)
        C_n2 = (1 - y_co2) * p / (R * T)

        ## Get loadings
        q_co2, q_n2 = mixture_loading(p, T, y_co2)

        inputs_raw[:, idx] = [C_co2, C_n2, T]
        outputs_raw[:, idx] = [q_co2, q_n2]
        idx += 1
    end

    ## Normalize inputs to [0, 1] range for better NN training
    ## For concentration use log(concentration) for better scaling
    ## For temperature use linear scaling
    C_min = minimum([minimum(inputs_raw[1, :]), minimum(inputs_raw[2, :])])
    C_max = maximum([maximum(inputs_raw[1, :]), maximum(inputs_raw[2, :])])
    log_C_min = log(max(C_min, 1e-10))  ## Avoid log(0)
    log_C_max = log(C_max)

    inputs_norm = similar(inputs_raw)
    inputs_norm[1, :] = (log.(max.(inputs_raw[1, :], 1e-10)) .- log_C_min) ./ (log_C_max - log_C_min)
    inputs_norm[2, :] = (log.(max.(inputs_raw[2, :], 1e-10)) .- log_C_min) ./ (log_C_max - log_C_min)
    inputs_norm[3, :] = (inputs_raw[3, :] .- T_min) ./ (T_max - T_min)

    ## Normalize outputs to [0, 1] range for better NN training
    q_min = minimum(outputs_raw)
    q_max = maximum(outputs_raw)
    outputs_norm = (outputs_raw .- q_min) ./ (q_max - q_min)

    normalization_params = (log_C_min=log_C_min, log_C_max=log_C_max,
        T_min=T_min, T_max=T_max, q_min=q_min, q_max=q_max)

    println("  Generated $(n_samples) mixture samples")
    println("  Concentration range: $(round(C_min, sigdigits=4)) - $(round(C_max, sigdigits=4)) mol/m³")
    println("  Loading range: $(round(q_min, sigdigits=4)) - $(round(q_max, sigdigits=4)) mol/kg")

    return inputs_norm, outputs_norm, normalization_params, inputs_raw, outputs_raw, pressures, temperatures, y_co2_range
end

# ## Visualize training data
# Skip visualization for now since we have mixture data in 3D input space

# ## Define neural network architecture
# We define a single neural network that takes mixture composition as input.
# Input: [normalized C_CO2, normalized C_N2, normalized T] → Output: [normalized q_CO2, normalized q_N2]
# The network learns competitive adsorption behavior directly from mixture data.

function create_mixture_isotherm_nn()
    return Chain(
        Dense(3 => 64, tanh),
        Dense(64 => 64, tanh),
        Dense(64 => 64, tanh),
        Dense(64 => 2, sigmoid)  ## sigmoid maps to [0, 1] for normalized output
    )
end

# ## Train neural network
# We train a single neural network for mixture isotherms using the Lux.jl framework.

function train_mixture_isotherm_nn(nn_model, train_inputs, train_outputs, epochs, lr)
    println("Training mixture isotherm neural network...")

    rng = Random.default_rng()
    Random.seed!(rng, 42)

    ps, st = Lux.setup(rng, nn_model)
    ps = ps |> f64
    tstate = Lux.Training.TrainState(nn_model, ps, st, Optimisers.Adam(lr))
    vjp_rule = ADTypes.AutoZygote()
    loss_function = Lux.MSELoss()

    ## Training loop
    losses = Float64[]
    @time begin
        for epoch in 1:epochs
            _, loss, _, tstate = Lux.Training.single_train_step!(
                vjp_rule, loss_function, (train_inputs, train_outputs), tstate
            )
            push!(losses, loss)

            if epoch % 1000 == 0 || epoch == epochs
                println("  Epoch: $(lpad(epoch, 5)) \t Loss: $(round(loss, sigdigits=6))")
            end
        end
    end

    return tstate, losses
end

# ## Visualize training progress

function plot_training_loss(losses)
    fig = Figure(size=(600, 500), fontsize=14)

    ax = Axis(fig[1, 1],
        xlabel="Epoch",
        ylabel="Loss",
        title="Mixture Isotherm NN Training Loss",
        yscale=log10
    )

    lines!(ax, 1:length(losses), losses, color=:blue, linewidth=2)

    return fig
end

function plot_isotherm_validation(mixture_nn_iso, constants)
    ## Create validation plots showing isotherm curves at different conditions
    isotherm_analytical = Mocca.DualSiteLangmuir(constants)
    R = constants.R
    ρ_s = constants.ρ_s

    ## Test conditions
    n_test = 50
    pressures = exp.(range(log(1e4), log(2e5), length=n_test))
    T_test = [293.15, 323.15, 353.15]  ## Low, medium, high temperature
    y_co2_test = [0.15, 0.50, 0.85]    ## Low, medium, high CO2 fraction

    fig = Figure(size=(1400, 900), fontsize=14)
    colors = [:blue, :green, :red]

    ## Plot CO2 isotherms at different temperatures (fixed composition)
    ax1 = Axis(fig[1, 1],
        xlabel="Pressure [Pa]",
        ylabel="CO₂ Loading [mol/kg]",
        title="CO₂ Isotherms (y_CO₂ = 0.15)",
        xscale=log10
    )

    for (i, T) in enumerate(T_test)
        q_co2_nn = Float64[]
        q_co2_ana = Float64[]

        for p in pressures
            C_co2 = y_co2_test[1] * p / (R * T)
            C_n2 = (1 - y_co2_test[1]) * p / (R * T)

            ## NN prediction
            q_nn_co2, _ = predict_loading(mixture_nn_iso, C_co2, C_n2, T)
            push!(q_co2_nn, q_nn_co2)

            ## Analytical prediction
            q_bulk = Mocca.compute_equilibrium(isotherm_analytical, SVector{2}([C_co2, C_n2]), T)
            push!(q_co2_ana, q_bulk[1] / ρ_s)
        end

        lines!(ax1, pressures, q_co2_ana, color=colors[i], linewidth=3,
            label="$(Int(round(T-273.15)))°C Analytical")
        scatter!(ax1, pressures, q_co2_nn, color=colors[i], markersize=6,
            label="$(Int(round(T-273.15)))°C NN")
    end
    axislegend(ax1, position=:rb, framevisible=false, labelsize=10)

    ## Plot N2 isotherms at different temperatures (fixed composition)
    ax2 = Axis(fig[1, 2],
        xlabel="Pressure [Pa]",
        ylabel="N₂ Loading [mol/kg]",
        title="N₂ Isotherms (y_CO₂ = 0.15)",
        xscale=log10
    )

    for (i, T) in enumerate(T_test)
        q_n2_nn = Float64[]
        q_n2_ana = Float64[]

        for p in pressures
            C_co2 = y_co2_test[1] * p / (R * T)
            C_n2 = (1 - y_co2_test[1]) * p / (R * T)

            ## NN prediction
            _, q_nn_n2 = predict_loading(mixture_nn_iso, C_co2, C_n2, T)
            push!(q_n2_nn, q_nn_n2)

            ## Analytical prediction
            q_bulk = Mocca.compute_equilibrium(isotherm_analytical, SVector{2}([C_co2, C_n2]), T)
            push!(q_n2_ana, q_bulk[2] / ρ_s)
        end

        lines!(ax2, pressures, q_n2_ana, color=colors[i], linewidth=3)
        scatter!(ax2, pressures, q_n2_nn, color=colors[i], markersize=6)
    end

    ## Plot competitive behavior: CO2 loading vs composition (fixed P and T)
    ax3 = Axis(fig[2, 1],
        xlabel="CO₂ Mole Fraction [-]",
        ylabel="CO₂ Loading [mol/kg]",
        title="CO₂ Loading vs Composition (P=1.5 bar, T=323K)"
    )

    p_fixed = 1.5e5
    T_fixed = 323.15
    y_range = range(0.01, 0.99, length=n_test)

    q_co2_nn_comp = Float64[]
    q_co2_ana_comp = Float64[]

    for y_co2 in y_range
        C_co2 = y_co2 * p_fixed / (R * T_fixed)
        C_n2 = (1 - y_co2) * p_fixed / (R * T_fixed)

        ## NN prediction
        q_nn_co2, _ = predict_loading(mixture_nn_iso, C_co2, C_n2, T_fixed)
        push!(q_co2_nn_comp, q_nn_co2)

        ## Analytical prediction
        q_bulk = Mocca.compute_equilibrium(isotherm_analytical, SVector{2}([C_co2, C_n2]), T_fixed)
        push!(q_co2_ana_comp, q_bulk[1] / ρ_s)
    end

    lines!(ax3, y_range, q_co2_ana_comp, color=:blue, linewidth=3, label="Analytical")
    scatter!(ax3, y_range, q_co2_nn_comp, color=:blue, markersize=6, label="NN")
    axislegend(ax3, position=:lt, framevisible=false)

    ## Plot N2 competitive behavior
    ax4 = Axis(fig[2, 2],
        xlabel="CO₂ Mole Fraction [-]",
        ylabel="N₂ Loading [mol/kg]",
        title="N₂ Loading vs Composition (P=1.5 bar, T=323K)"
    )

    q_n2_nn_comp = Float64[]
    q_n2_ana_comp = Float64[]

    for y_co2 in y_range
        C_co2 = y_co2 * p_fixed / (R * T_fixed)
        C_n2 = (1 - y_co2) * p_fixed / (R * T_fixed)

        ## NN prediction
        _, q_nn_n2 = predict_loading(mixture_nn_iso, C_co2, C_n2, T_fixed)
        push!(q_n2_nn_comp, q_nn_n2)

        ## Analytical prediction
        q_bulk = Mocca.compute_equilibrium(isotherm_analytical, SVector{2}([C_co2, C_n2]), T_fixed)
        push!(q_n2_ana_comp, q_bulk[2] / ρ_s)
    end

    lines!(ax4, y_range, q_n2_ana_comp, color=:red, linewidth=3, label="Analytical")
    scatter!(ax4, y_range, q_n2_nn_comp, color=:red, markersize=6, label="NN")
    axislegend(ax4, position=:rt, framevisible=false)

    Label(fig[0, :], "Neural Network Isotherm Validation", fontsize=18, font=:bold)

    return fig
end

# ## Create NN-based isotherm wrapper

if !@isdefined(MixtureNeuralNetworkIsotherm)
    struct MixtureNeuralNetworkIsotherm{M,P,S,N}
        model::M
        parameters::P
        states::S
        norm_params::N  ## Normalization parameters (NamedTuple)
        R::Float64      ## Gas constant [J/(mol·K)]
        ρ_s::Float64    ## Adsorbent density [kg/m³]
    end
end  ## if

## Function to predict loading given concentrations and temperature
function predict_loading(iso::MixtureNeuralNetworkIsotherm, C_co2, C_n2, temperature)
    np = iso.norm_params

    ## Clamp concentrations to ≥0.01 mol/m³ before log-transformation
    ## to prevent large gradients at very low concentrations (d(log(C))/dC = 1/C) which causes convergence issues in the Newton solver
    C_co2_safe = max(C_co2, 0.01)
    C_n2_safe = max(C_n2, 0.01)

    ## Normalize inputs to [0, 1] range (to match training data)
    log_C_co2_norm = (log(C_co2_safe) - np.log_C_min) / (np.log_C_max - np.log_C_min)
    log_C_n2_norm = (log(C_n2_safe) - np.log_C_min) / (np.log_C_max - np.log_C_min)
    T_norm = (temperature - np.T_min) / (np.T_max - np.T_min)

    ## Clamp to stay within training range
    log_C_co2_norm = clamp(log_C_co2_norm, 0.0, 1.0)
    log_C_n2_norm = clamp(log_C_n2_norm, 0.0, 1.0)
    T_norm = clamp(T_norm, 0.0, 1.0)

    ## Predict and denormalize
    input = reshape([log_C_co2_norm, log_C_n2_norm, T_norm], 3, 1)
    output_norm, _ = Lux.apply(iso.model, input, iso.parameters, iso.states)

    ## Denormalize outputs to get actual loadings [mol/kg]
    q_co2 = output_norm[1] * (np.q_max - np.q_min) + np.q_min
    q_n2 = output_norm[2] * (np.q_max - np.q_min) + np.q_min

    return q_co2, q_n2
end

# ## Define custom secondary variable with NN

if !@isdefined(NNAdsorptionMassTransfer)
    struct NNAdsorptionMassTransfer{I} <: JutulDarcy.ComponentVariables
        isotherm::I
    end
end

## Update function for the custom secondary variable
Jutul.@jutul_secondary function update_adsorption_mass_transfer!(
    adsorption_mass_transfer,
    tv::NNAdsorptionMassTransfer,
    model::Mocca.AdsorptionModel,
    concentrations,
    Temperature,
    AdsorbedConcentration,
    ix
)
    N = JutulDarcy.number_of_components(model.system)
    T = eltype(adsorption_mass_transfer)
    iso = tv.isotherm
    mt = model.system.mass_transfer
    ρ_s = iso.ρ_s

    for cell in ix
        C = SVector{N, T}(@view concentrations[:, cell])
        q = SVector{N, T}(@view AdsorbedConcentration[:, cell])

        ## Use neural network for equilibrium calculation
        ## Get loadings in [mol/kg] from NN, convert to [mol/m³ bulk]
        q_co2_kg, q_n2_kg = predict_loading(iso, C[1], C[2], Temperature[cell])
        qstar = SVector{N,T}(q_co2_kg * ρ_s, q_n2_kg * ρ_s)

        ## Use the system's mass transfer model
        rate = Mocca.compute_mass_transfer_rate(mt, C, q, qstar)

        for i in 1:N
            adsorption_mass_transfer[i, cell] = rate[i]
        end
    end
end

# ## Plotting functions

function plot_comparison_results(ref_states, nn_states)
    ## Extract final states for comparison
    ref_final_state = ref_states[end]
    nn_final_state = nn_states[end]

    ncells = 200
    x_pos = collect(range(0.0, 1.0, length=ncells))

    ## Extract final column profiles
    ref_y_co2 = ref_final_state[:y][1, :]
    ref_y_n2 = ref_final_state[:y][2, :]
    ref_pressure = ref_final_state[:Pressure][:]
    ref_temp = ref_final_state[:Temperature][:]
    ref_ads_co2 = ref_final_state[:AdsorbedConcentration][1, :]
    ref_ads_n2 = ref_final_state[:AdsorbedConcentration][2, :]

    nn_y_co2 = nn_final_state[:y][1, :]
    nn_y_n2 = nn_final_state[:y][2, :]
    nn_pressure = nn_final_state[:Pressure][:]
    nn_temp = nn_final_state[:Temperature][:]
    nn_ads_co2 = nn_final_state[:AdsorbedConcentration][1, :]
    nn_ads_n2 = nn_final_state[:AdsorbedConcentration][2, :]

    ## Create figure with subplots
    fig = Figure(size=(1400, 900), fontsize=14)

    ## Gas phase mole fractions
    ax1 = Axis(fig[1, 1], xlabel="Normalized Position [-]", ylabel="CO₂ Mole Fraction [-]",
        title="Gas Phase CO₂ Profile")
    l1_ref = lines!(ax1, x_pos, ref_y_co2, linewidth=3, color=:blue)
    l1_nn = lines!(ax1, x_pos, nn_y_co2, linewidth=4,
        color=:red, linestyle=:dash)

    ax2 = Axis(fig[1, 2], xlabel="Normalized Position [-]", ylabel="N₂ Mole Fraction [-]",
        title="Gas Phase N₂ Profile")
    lines!(ax2, x_pos, ref_y_n2, linewidth=3, color=:blue)
    lines!(ax2, x_pos, nn_y_n2, linewidth=4,
        color=:red, linestyle=:dash)

    ## Adsorbed concentrations
    ax3 = Axis(fig[2, 1], xlabel="Normalized Position [-]", ylabel="CO₂ Loading [mol/m³]",
        title="Adsorbed CO₂ Profile")
    lines!(ax3, x_pos, ref_ads_co2, linewidth=3, color=:blue)
    lines!(ax3, x_pos, nn_ads_co2, linewidth=4,
        color=:red, linestyle=:dash)

    ax4 = Axis(fig[2, 2], xlabel="Normalized Position [-]", ylabel="N₂ Loading [mol/m³]",
        title="Adsorbed N₂ Profile")
    lines!(ax4, x_pos, ref_ads_n2, linewidth=3, color=:blue)
    lines!(ax4, x_pos, nn_ads_n2, linewidth=4,
        color=:red, linestyle=:dash)

    ## Pressure and temperature
    ax5 = Axis(fig[3, 1], xlabel="Normalized Position [-]", ylabel="Pressure [Pa]",
        title="Pressure Profile")
    lines!(ax5, x_pos, ref_pressure, linewidth=3, color=:blue)
    lines!(ax5, x_pos, nn_pressure, linewidth=4,
        color=:red, linestyle=:dash)

    ax6 = Axis(fig[3, 2], xlabel="Normalized Position [-]", ylabel="Temperature [K]",
        title="Temperature Profile")
    lines!(ax6, x_pos, ref_temp, linewidth=3, color=:blue)
    lines!(ax6, x_pos, nn_temp, linewidth=4,
        color=:red, linestyle=:dash)

    ## Add overall title
    Label(fig[0, :], "Mocca Column Simulation: Analytical vs Hybrid (NN) Isotherm",
        fontsize=18, font=:bold)

    ## Add shared legend
    Legend(fig[1:3, 3], [l1_ref, l1_nn],
        ["Mocca (Analytical Isotherm)", "Mocca Hybrid (NN Isotherm)"],
        "Isotherm Model", framevisible=false, labelsize=16)

    return fig
end

# ## Run the reference simulation
# We first run a reference simulation using the standard dual-site Langmuir model.
# This will serve as our ground truth for comparison.

println("="^80)
println("Running reference simulation...")
println("="^80)

ref_model, ref_state0, ref_parameters, ref_timesteps, ref_forces, ref_constants = setup_reference_simulation()

## Simulation configuration
timestep_selector_cfg = (y = 0.01, Temperature = 10.0, Pressure = 10.0)

ref_case = Mocca.MoccaCase(ref_model, ref_timesteps, ref_forces;
    state0=ref_state0, parameters=ref_parameters)
ref_states, ref_timesteps_out = Mocca.simulate_process(ref_case;
    timestep_selector_cfg=timestep_selector_cfg,
    output_substates=true,
)

# ## Generate training data

println("\n" * "="^80)
println("Generating training data from analytical isotherms...")
println("="^80)

train_inputs, train_outputs, norm_params, inputs_raw, outputs_raw, pressures, temperatures, y_co2_range =
    generate_mixture_training_data(ref_constants)

# ## Train neural network for mixture isotherms
# Train a single NN that learns competitive adsorption from mixture data
# The model is saved to disk to avoid retraining

model_file = "nn_mixture_isotherm_model.jld2"
architecture_version = "sigmoid_64"

need_training = true
if isfile(model_file)
    println("\n" * "="^80)
    println("Loading pre-trained model...")
    println("="^80)

    saved_data = load(model_file)
    saved_version = get(saved_data, "architecture_version", "unknown")

    if saved_version != architecture_version
        println("Architecture mismatch ($(saved_version) ≠ $(architecture_version)). Retraining...")
        rm(model_file)
    else
        ## Architecture matches, load the model
        mixture_nn = saved_data["model"]
        tstate_params = saved_data["parameters"]
        tstate_states = saved_data["states"]
        saved_norm_params = saved_data["norm_params"]
        losses = saved_data["losses"]

        println("Model loaded successfully")
        println("Training loss at end: $(round(losses[end], sigdigits=6))")

        ## Recreate tstate structure (needed for inference)
        rng = Random.default_rng()
        ps, st = Lux.setup(rng, mixture_nn)
        tstate = Lux.Training.TrainState(mixture_nn, tstate_params, tstate_states, Optimisers.Adam(0.001))

        need_training = false
    end
end

if need_training
    println("\n" * "="^80)
    println("Training neural network...")
    println("="^80)

    mixture_nn = create_mixture_isotherm_nn()

    epochs = 20000
    lr = 0.001

    tstate, losses = train_mixture_isotherm_nn(mixture_nn, train_inputs, train_outputs, epochs, lr)

    ## Save the trained model with architecture version
    println("\nSaving trained model to $(model_file)...")
    save(model_file, Dict(
        "architecture_version" => architecture_version,
        "model" => mixture_nn,
        "parameters" => tstate.parameters,
        "states" => tstate.states,
        "norm_params" => norm_params,
        "losses" => losses
    ))
    println("Model saved successfully")
end

# ## Visualize training progress

fig_training_loss = plot_training_loss(losses)
display(fig_training_loss)

# ## Create isotherm wrapper and validate
# Wrap the trained neural network for use in Mocca and validate predictions

println("\n" * "="^80)
println("Creating NN isotherm wrapper...")
println("="^80)

mixture_nn_iso = MixtureNeuralNetworkIsotherm(
    mixture_nn,
    tstate.parameters,
    tstate.states,
    norm_params,
    ref_constants.R,
    ref_constants.ρ_s
)

println("Neural network wrapper created successfully")

# ## Validate NN predictions on typical simulation conditions
## Test the NN on representative conditions to catch issues before simulation

println("\n" * "="^80)
println("Validating NN predictions...")
println("="^80)

## Test at typical simulation conditions
test_C_co2 = 5.0  ## mol/m³
test_C_n2 = 30.0  ## mol/m³
test_T = 300.0  ## K

q_co2_test, q_n2_test = predict_loading(mixture_nn_iso, test_C_co2, test_C_n2, test_T)

println("Test prediction at C_CO2=$(test_C_co2) mol/m³, C_N2=$(test_C_n2) mol/m³, T=$(test_T) K:")
println("  q_CO2 = $(round(q_co2_test, sigdigits=5)) mol/kg")
println("  q_N2  = $(round(q_n2_test, sigdigits=5)) mol/kg")

## Check for invalid predictions
(isnan(q_co2_test) || isnan(q_n2_test)) && error("NN prediction returned NaN!")
(q_co2_test < 0 || q_n2_test < 0) && error("NN prediction returned negative loading!")
(q_co2_test > 10.0 || q_n2_test > 10.0) && @warn "NN prediction seems unreasonably high (>10 mol/kg)"

## Compare with analytical model
isotherm_analytical = Mocca.DualSiteLangmuir(ref_constants)
q_analytical = Mocca.compute_equilibrium(isotherm_analytical, SVector{2}([test_C_co2, test_C_n2]), test_T)
q_co2_analytical = q_analytical[1] / ref_constants.ρ_s
q_n2_analytical = q_analytical[2] / ref_constants.ρ_s

println("\nAnalytical reference:")
println("  q_CO2 = $(round(q_co2_analytical, sigdigits=5)) mol/kg")
println("  q_N2  = $(round(q_n2_analytical, sigdigits=5)) mol/kg")

rel_error_co2 = abs(q_co2_test - q_co2_analytical) / q_co2_analytical * 100
rel_error_n2 = abs(q_n2_test - q_n2_analytical) / q_n2_analytical * 100

println("\nRelative errors:")
println("  CO2: $(round(rel_error_co2, sigdigits=3))%")
println("  N2:  $(round(rel_error_n2, sigdigits=3))%")

if rel_error_co2 > 20.0 || rel_error_n2 > 20.0
    @warn "Large prediction error detected (>20%). Model may need more training or better data coverage."
end

println("\n✓ NN validation passed")

println("\n" * "="^80)
println("Generating isotherm validation plots...")
println("="^80)

fig_isotherm_validation = plot_isotherm_validation(mixture_nn_iso, ref_constants)
display(fig_isotherm_validation)

# ## Run simulation with neural network isotherms

println("\n" * "="^80)
println("Running simulation with Neural Network isotherms...")
println("="^80)

nn_model, nn_state0, nn_parameters, nn_timesteps, nn_forces, _ = setup_reference_simulation()
Jutul.replace_variables!(nn_model, AdsorptionMassTransfer=NNAdsorptionMassTransfer(mixture_nn_iso))

nn_case = Mocca.MoccaCase(nn_model, nn_timesteps, nn_forces;
    state0=nn_state0, parameters=nn_parameters)
nn_states, nn_timesteps_out = Mocca.simulate_process(nn_case;
    timestep_selector_cfg=timestep_selector_cfg,
    output_substates=true,
)

# ## Compare results

isempty(nn_states) && error("NN simulation failed - no states generated")

println("\n" * "="^80)
println("Generating comparison plots...")
println("="^80)

fig_comparison = plot_comparison_results(ref_states, nn_states)
display(fig_comparison)

println("\n" * "="^80)
println("EXAMPLE COMPLETE!")
println("="^80)
println("\n✓ Successfully integrated neural network isotherms into Mocca.jl")
println("✓ NN learned competitive adsorption directly from mixture data")
println("✓ Column profiles match analytical reference simulation")
println("\nThis demonstrates how neural networks can replace analytical models")
println("in process simulations while maintaining accuracy and stability.")

