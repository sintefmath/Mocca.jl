# AdjointsDI: `deps = :case` gives wrong gradients when the first report step is split

**Package:** Jutul 0.4.25 and 0.4.29, `src/ad/AdjointsDI/adjoints.jl`
**Reproduction:** [`jutul_adjointsdi_split_step.jl`](jutul_adjointsdi_split_step.jl) (Jutul only, about a minute)

## Summary

`solve_adjoint_generic` with `deps = :case` gives wrong gradients when the simulator splits report step 1 into several sub-steps. This is also the default path for `DictOptimization.parameters_gradient`. Splitting any later report step is handled correctly, and `deps = :parameters` is correct in all cases.

## Cause

`evaluate_residual_and_jacobian_for_state_pair` contains:

```julia
if step_info[:step] == 1
    state0 = case.state0
end
```

`step_info[:step]` is the report step, not the sub-step. Every sub-step of report step 1 therefore has its residual linearised from the initial state, rather than from the previous sub-step that `adjoint_step_state_triplet` passed in. The override is only needed for the first sub-step: there, `case.state0` carries the dependence on `x` through the setup function. The fix is:

```julia
if step_info[:substep_global] == 1
    state0 = case.state0
end
```

## Reproduction

The script uses the time-dependent `VariablePoissonSystem` case from `test/adjoints/basic_adjoint.jl`, with `x = [dx, dy, U0, k_val, srcval]`. It sets `max_timestep` below the report step length so that the simulator splits steps, and runs with `output_substates = true`. It compares the adjoint gradient against central finite differences for three scenarios:

- **A:** one report step, not split.
- **B:** one report step, split into four.
- **C:** two report steps, with only the second one split.

It then applies the fix above in place and reruns the scenarios.

The objective is `dt * U[end]^2`. The objective in Jutul's own test, `U[end] - U[1]`, hides the bug: it doesn't change under a uniform shift of `U0`, and in this model `state0` only enters the gradient through `U0`.

```
Jutul 0.4.29, before patch:
A: one step, not split             sub-steps per report step [1]         :case rel. error 4.2e-09   :parameters rel. error 1.3e-08
B: first report step split         sub-steps per report step [4]         :case rel. error 1.5e+00   :parameters rel. error 1.2e-09
C: only second report step split   sub-steps per report step [1, 4]      :case rel. error 5.8e-10   :parameters rel. error 1.8e-09

After patch (step_info[:step] == 1  ->  step_info[:substep_global] == 1):
A: one step, not split             sub-steps per report step [1]         :case rel. error 4.2e-09   :parameters rel. error 1.3e-08
B: first report step split         sub-steps per report step [4]         :case rel. error 5.8e-10   :parameters rel. error 1.2e-09
C: only second report step split   sub-steps per report step [1, 4]      :case rel. error 5.8e-10   :parameters rel. error 1.8e-09
```

In scenario B, only the `U0` entry is wrong: 666.4 against the correct 266.4. The other entries don't interact with `state0` in this model, because its accumulation term has no cell volume. In a model whose accumulation term depends on the parameters, every such parameter is affected.

## Where it showed up

The problem appeared in Mocca.jl, a Jutul-based adsorption process simulator. The gradient of a vacuum-product objective with respect to column length came out 1.2% off with `DictOptimization.parameters_gradient`, using the default `deps = :case`. That run splits only report step 1 (21 sub-steps for 20 report steps). Applying the same one-line change brings the error to 5e-11.

## Suggested test

Add a split first step to `test_for_timesteps` in `test/adjoints/basic_adjoint.jl`:
- Pass `max_timestep` through to `simulate`, both in the forward run and in `num_grad_generic`.
- Use `output_substates = true`.
- Use an objective that depends on `U0`, for example `U[end]^2`.
