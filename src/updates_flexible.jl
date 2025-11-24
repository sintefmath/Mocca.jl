# Extension of Mocca update functions to support flexible isotherm models
# The FlexibleAdsorptionSystem delegates isotherm calculations to pluggable isotherm models

# Override compute_equilibrium for FlexibleAdsorptionSystem to use the pluggable isotherm model
function compute_equilibrium(sys::FlexibleAdsorptionSystem, concentration, temperature)
    return compute_equilibrium(sys.isotherm_model, concentration, temperature)
end

# Override compute_ki for FlexibleAdsorptionSystem to use the pluggable isotherm model
function compute_ki(sys::FlexibleAdsorptionSystem, concentration, qstar)
    return compute_ki(sys.isotherm_model, concentration, qstar)
end

# Note: The existing update_adsorption_mass_transfer function in updates.jl
# automatically works with FlexibleAdsorptionSystem through multiple dispatch,
# as it calls compute_equilibrium and compute_ki which are now overridden above.
