# Minimal reproduction: Jutul.AdjointsDI.solve_adjoint_generic with deps = :case
# gives wrong gradients when the first report step is split into sub-steps.
# deps = :parameters is correct in the same situation.
#
# Uses only Jutul (the time-dependent VariablePoissonSystem from Jutul's own
# adjoint tests). Run with:  julia --project=<env with Jutul> jutul_adjointsdi_split_step.jl
using Jutul, Printf
import Jutul.AdjointsDI: solve_adjoint_generic

# Same test case as Jutul's test/adjoints/basic_adjoint.jl.
# x = [dx, dy, U0, k_val, srcval]
function setup_poisson_case(x, step_info = missing; dt = [1.0], dim = (3, 1))
    dx, dy, U0, k_val, srcval = x
    sys = VariablePoissonSystem(time_dependent = true)
    g = CartesianMesh(dim, (dx, dy))
    D = DiscretizedDomain(g, (poisson = Jutul.PoissonDiscretization(g),))
    model = SimulationModel(D, sys)
    state0 = setup_state(model, Dict(:U => U0))
    parameters = setup_parameters(model, K = compute_face_trans(g, k_val))
    nc = number_of_cells(g)
    forces = setup_forces(model, sources = [PoissonSource(1, srcval), PoissonSource(nc, -srcval)])
    return JutulCase(model, dt, forces; parameters = parameters, state0 = state0)
end

G(model, state, dt, step_info, forces) = dt * state[:U][end]^2

function objective(x; dt, max_timestep)
    case = setup_poisson_case(x; dt = dt)
    res = simulate(case, info_level = -1, max_timestep = max_timestep, output_substates = true)
    states, dts, = Jutul.expand_to_ministeps(res)
    return sum(G(case.model, s, h, nothing, nothing) for (s, h) in zip(states, dts))
end

function fd_gradient(x0; kwarg...)
    h = 1e-6
    return map(eachindex(x0)) do i
        xp = copy(x0); xp[i] += h
        xm = copy(x0); xm[i] -= h
        (objective(xp; kwarg...) - objective(xm; kwarg...)) / (2h)
    end
end

function adjoint_gradient(x0; dt, max_timestep, deps)
    case = setup_poisson_case(x0; dt = dt)
    res = simulate(case, info_level = -1, max_timestep = max_timestep, output_substates = true)
    F = (x, step_info) -> setup_poisson_case(x, step_info; dt = dt)
    return solve_adjoint_generic(x0, F, res.states, res.reports, G;
        state0 = case.state0, forces = case.forces, deps = deps, deps_ad = :jutul, info_level = -1)
end

function run_scenario(name; dt, max_timestep)
    x0 = [1.0, 1.0, 1.0, 1.0, 1.0]
    case = setup_poisson_case(x0; dt = dt)
    res = simulate(case, info_level = -1, max_timestep = max_timestep, output_substates = true)
    _, _, report_ix = Jutul.expand_to_ministeps(res)
    nsub = [count(==(i), report_ix) for i in eachindex(dt)]
    g_fd = fd_gradient(x0; dt = dt, max_timestep = max_timestep)
    g_case = adjoint_gradient(x0; dt = dt, max_timestep = max_timestep, deps = :case)
    g_prm = adjoint_gradient(x0; dt = dt, max_timestep = max_timestep, deps = :parameters)
    # deps = :parameters only covers dx, dy and k_val (entries 1, 2, 4).
    ix = [1, 2, 4]
    relerr(a, b) = maximum(abs.(a .- b)) / maximum(abs.(b))
    @printf("%-34s sub-steps per report step %-10s  :case rel. error %.1e   :parameters rel. error %.1e\n",
        name, string(nsub), relerr(g_case, g_fd), relerr(g_prm[ix], g_fd[ix]))
    return (fd = g_fd, case = g_case, parameters = g_prm)
end

function run_all()
    run_scenario("A: one step, not split"; dt = [100.0], max_timestep = Inf)
    run_scenario("B: first report step split"; dt = [100.0], max_timestep = 25.0)
    run_scenario("C: only second report step split"; dt = [10.0, 100.0], max_timestep = 25.0)
end

println("Jutul ", pkgversion(Jutul), ", before patch:")
run_all()

# Cause, src/ad/AdjointsDI/adjoints.jl, evaluate_residual_and_jacobian_for_state_pair:
#
#     if step_info[:step] == 1
#         state0 = case.state0
#     end
#
# step_info[:step] is the report step, so every sub-step of report step 1 is
# linearised from the initial state instead of the previous sub-step. The state
# triplet already supplies the right state0 (packed_steps.state0 for the first
# sub-step), so the override is only needed for the first sub-step, where it
# makes state0 depend on x:
#
#     if step_info[:substep_global] == 1
#
# Apply that change in place and rerun. run_all is called with invokelatest so
# that it sees the redefined method.
if get(ENV, "APPLY_PATCH", "1") == "1"
    file = joinpath(pkgdir(Jutul), "src", "ad", "AdjointsDI", "adjoints.jl")
    src = read(file, String)
    old = "if step_info[:step] == 1\n        state0 = case.state0"
    occursin(old, src) || error("Could not find the line to patch in $file")
    patched = replace(src, old => "if step_info[:substep_global] == 1\n        state0 = case.state0")
    Base.include_string(Jutul.AdjointsDI, patched, file)
    println("\nAfter patch (step_info[:step] == 1  ->  step_info[:substep_global] == 1):")
    Base.invokelatest(run_all)
end
