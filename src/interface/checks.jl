# Checks of the gradients of a process against finite differences, at two step
# sizes so that the error of the steps can be told apart from a gradient that
# is wrong.

# Start the next steady state of the process from result `r`
function _warm_from!(p::CyclicProcess, r)
    p.warm.state = r.state0
    p.warm.jacobian = r.jacobian
    return p
end
function _warm_from!(p::SteadyProcess, r)
    p.warm.state = restart_state(p.model, r.states[end])
    return p
end

# The state the process settles to, converged tightly enough for finite
# differences, from `start` (a result of an earlier solve) if given
function _tight_steady_state(p::Process, start; tol)
    isnothing(start) || _warm_from!(p, start)
    return p isa CyclicProcess ? cyclic_steady_state(p; tol = tol, progress = false) :
        steady_state(p; tol = tol, max_spans = 100, progress = false)
end

"""
    check_gradients(p; variables = all, measures = all, steps = (1e-3, 1e-4), rtol = 1e-3,
        tol = 1e-10, io = stdout) → NamedTuple

Compare the gradients of the performance of process `p` from
[`performance_gradient`](@ref) with central differences of its performance at
the state it settles to, converged to `tol`, for each of `variables` and
`measures`. Each variable is changed by each of the relative `steps` (an
absolute step for a variable at zero), each difference starting from the
same steady state.

Each pair of variable and measure is
- `:pass` if the gradient agrees with one of the differences to `rtol`;
- `:inconclusive` if not, but the differences at the two steps disagree with
  each other as much: the performance has a kink or a jump there, such as a
  pump that switches on, or the steady states are not converged enough for
  the steps;
- `:fail` otherwise: the differences agree with each other but not with the
  gradient, so the gradient is wrong. A model that breaks one of the rules of
  the extending guide (see the docs, *Rules for differentiable models*) ends
  here.

Agreement is relative to the larger of the two, or, for a derivative near
zero, to the change of the measure over the variable's value. Prints a table
to `io` (`nothing` for none) and returns the `rows` (`variable`, `measure`,
`gradient`, `differences`, `error`, `status`), whether every row `passed`
(none failed) and the `metrics` at the steady state. Leaves the process as it
was, with its last steady state.
"""
function check_gradients(p::Process; variables = variable_names(p), measures = keys(performance_measures(p)),
        steps = (1e-3, 1e-4), rtol = 1e-3, tol = 1e-10, io = stdout)
    names = collect(Symbol, variables)
    ms = collect(Symbol, measures)
    isempty(names) && input_error("check_gradients", "variables", "give at least one design variable")
    x0 = variable_values(p, names)
    base = _tight_steady_state(p, nothing; tol = tol)
    base.converged === true || @warn "The process did not settle to $tol; the checks may be inconclusive"
    g = performance_gradient(p; variables = names, measures = ms, steady_state = base, progress = false)
    rows = Any[]
    for (i, n) in enumerate(names)
        differences = map(steps) do s
            h = x0[i] == 0 ? s : s * abs(x0[i])
            perf(x) = (p[n] = x; r = _tight_steady_state(p, base; tol = tol); m = performance(r); p[n] = x0[i]; (metrics = m, settled = r.converged === true))
            up, down = perf(x0[i] + h), perf(x0[i] - h)
            (up.settled && down.settled) || @warn "The process did not settle at $n ± $h"
            map(m -> (up.metrics[m] - down.metrics[m]) / 2h, NamedTuple{Tuple(ms)}(Tuple(ms)))
        end
        for m in ms
            a = g.gradients[m][n]
            d = [δ[m] for δ in differences]
            # A derivative near zero is measured against the change of the
            # measure over the variable's value
            err, status = _check_status(a, d, 1e-6 * abs(g.metrics[m]) / max(abs(x0[i]), eps()), rtol)
            push!(rows, (variable = n, measure = m, gradient = a, differences = Tuple(d), error = err, status = status))
        end
    end
    _warm_from!(p, base)
    isnothing(io) || _print_checks(io, p, rows, steps)
    return (rows = rows, passed = all(r -> r.status != :fail, rows), metrics = g.metrics)
end

# The error of `gradient` against the closest of the `differences`, relative
# to the larger of the two or `floor`, and whether it passes, fails, or the
# differences disagree with each other more than with the gradient
function _check_status(gradient, differences, floor, rtol)
    rel(u, v) = abs(u - v) / max(abs(u), abs(v), floor)
    err = minimum(rel(gradient, δ) for δ in differences)
    spread = maximum(rel(a, b) for a in differences, b in differences)
    return (err, err <= rtol ? :pass : spread > err ? :inconclusive : :fail)
end

function _print_checks(io, p, rows, steps)
    println(io, "Gradient checks", isempty(p.name) ? "" : " of $(p.name)", ": central differences at relative steps ", join(steps, " and "))
    println(io, rpad("variable", 18), rpad("measure", 26), lpad("gradient", 13), join(lpad("step $s", 13) for s in steps), lpad("error", 10), "  status")
    for r in rows
        println(io, rpad(r.variable, 18), rpad(r.measure, 26), lpad(_sci(r.gradient), 13), join(lpad(_sci(d), 13) for d in r.differences),
            lpad(_sci(r.error), 10), "  ", r.status)
    end
    counts = [count(r -> r.status == s, rows) for s in (:pass, :inconclusive, :fail)]
    println(io, "$(counts[1]) passed, $(counts[2]) inconclusive, $(counts[3]) failed")
end
