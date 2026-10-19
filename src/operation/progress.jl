# ------------------------------------------------------------------------------
# Progress of the cyclic steady state solvers
# ------------------------------------------------------------------------------

"""
    CycleProgress(io = stdout; live = io isa Base.TTY)

Printout of the progress of the cyclic steady state solvers, given to them as
`progress`. Each solver is a phase with a tolerance of its own (such as
Anderson acceleration to `1e-3`, then Newton's method to the final tolerance),
and each phase has a line with a bar. The line shows the cycle count, the
measure of change against the phase's tolerance, a bar of how far the phase has
come (on a log scale, from its first cycle to its tolerance, by the smallest
change so far, as the change of an accelerated cycle goes up and down), an
estimate of the cycles left from the recent rate of decrease of that smallest
change, and the time taken. When a phase ends, its line says whether it
converged, or for a phase that hands over to the next, that it is done.

With `live`, as in a terminal, the lines of the phases planned with
`plan_phases!` are drawn together and rewritten in place, the later
ones waiting until their turn. Otherwise each cycle gets a line of its own.
"""
mutable struct CycleProgress
    io::IO
    live::Bool
    started::Float64
    # The phases planned, as (name, tol, noun), and the line of each
    plan::Vector{Tuple{String, Float64, String}}
    lines::Vector{String}
    # The current phase: its place in the plan, and its settings
    current::Int
    noun::String
    measure::String
    tol::Float64
    estimate::Bool
    phase_started::Float64
    values::Vector{Float64}
    # Time taken by each cycle of the phase
    durations::Vector{Float64}
    # Cycles simulated in every phase
    cycles::Int
    # Lines of the live block on screen, and when it was last drawn
    drawn::Int
    last_status::Float64
end

CycleProgress(io::IO = stdout; live = io isa Base.TTY) =
    CycleProgress(io, live, time(), Tuple{String, Float64, String}[], String[], 0, "cycle", "change", NaN, true,
        time(), Float64[], Float64[], 0, 0, 0.0)

"""
    cycle_progress(progress; info_level = -1) → CycleProgress or nothing

The printout for the `progress` keyword of the cyclic steady state solvers:
`false` for none, `true` for standard output, an `IO` to print to, or a
[`CycleProgress`](@ref) to continue. Lines are not rewritten in place while
Jutul prints its own (`info_level ≥ 0`).
"""
cycle_progress(progress::CycleProgress; info_level = -1) = progress
cycle_progress(progress::Bool; info_level = -1) = progress ? cycle_progress(stdout; info_level = info_level) : nothing
cycle_progress(progress::IO; info_level = -1) = CycleProgress(progress; live = progress isa Base.TTY && info_level < 0)
cycle_progress(::Nothing; info_level = -1) = nothing

# Every report is a no-op without a printout
plan_phases!(::Nothing, args...) = nothing
start_phase!(::Nothing, args...; kwarg...) = nothing
report_cycle!(::Nothing, args...; kwarg...) = nothing
report_status!(::Nothing, args...; kwarg...) = nothing
finish_phase!(::Nothing, args...; kwarg...) = nothing
report_note!(::Nothing, args...) = nothing
count_cycle!(::Nothing) = nothing
finish_progress!(::Nothing, args...) = nothing

"""
    plan_phases!(pr, phases)

Plan the phases to come, as `(name, tol, noun)`, each with its own line and
bar, so that a live printout shows the later ones waiting. A phase started
without a plan is planned on its own.
"""
function plan_phases!(pr::CycleProgress, phases)
    _close_block!(pr)
    pr.plan = [(String(n), Float64(t), String(u)) for (n, t, u) in phases]
    pr.lines = [_waiting_line(pr, k) for k in eachindex(pr.plan)]
    pr.current = 0
    pr.live && _draw_block!(pr)
    return pr
end

"""
    start_phase!(pr, name; tol, noun = "cycle", measure = "change", estimate = true)

Start a phase `name` that runs until its `measure` is below `tol`, counting
in `noun`s, with an estimate of the number left if `estimate`. A phase with
`tol = NaN` has no target, such as a set number of cycles. It takes the
next line of the plan if that is `name`, and otherwise plans it on its own.
"""
function start_phase!(pr::CycleProgress, name; tol, noun = "cycle", measure = "change", estimate = true)
    k = pr.current + 1
    if k > length(pr.plan) || pr.plan[k][1] != name
        plan_phases!(pr, [(name, tol, noun)])
        k = 1
    end
    pr.current = k
    pr.noun, pr.measure, pr.tol, pr.estimate = noun, measure, tol, estimate
    pr.phase_started = time()
    empty!(pr.values)
    empty!(pr.durations)
    return pr
end

"""
    report_cycle!(pr, n, value; left = nothing, fraction = nothing, note = "")

Report cycle (or iteration) `n` of the current phase, with its measure
`value`. `left` and `fraction` replace the estimate of the cycles left and of
how far the phase has come, and `note` is added to the line.
"""
function report_cycle!(pr::CycleProgress, n, value; left = nothing, fraction = nothing, note = "")
    push!(pr.values, value)
    push!(pr.durations, time() - pr.phase_started - sum(pr.durations))
    pr.cycles += 1
    left = isnothing(left) && pr.estimate ? _cycles_to_reach(pr.values, pr.tol) : left
    frac = something(fraction, _log_fraction(pr.values, pr.tol))
    target = isnan(pr.tol) ? "" : " → $(_sci(pr.tol))"
    parts = [_phase_columns(pr, pr.current, "$(pr.noun) $(lpad(n, 3))", "$(pr.measure) $(_sci(value))" * target, frac)]
    isempty(note) || push!(parts, note)
    if !isnothing(left) && left > 0
        # The last few cycles, as the first ones include compilation
        recent = pr.durations[max(end - 4, 1):end]
        per_cycle = sum(recent) / length(recent)
        push!(parts, "~$left $(_plural(pr.noun, left)), ~$(_duration(left * per_cycle)) left")
    end
    push!(parts, _duration(time() - pr.phase_started))
    _update_line!(pr, join(parts, "  "))
    return pr
end

"""
    report_status!(pr, message; persist = false)

Show what the current phase is doing between cycles, such as building a
Jacobian. Unless `persist`, it is only shown on a live printout, and at most
four times a second.
"""
function report_status!(pr::CycleProgress, message; persist = false)
    (persist || (pr.live && time() - pr.last_status > 0.25)) || return pr
    if pr.live
        # In place of the end of the phase's line
        line = pr.lines[pr.current]
        _update_line!(pr, first(line, min(length(line), _status_column())) * message)
    else
        println(pr.io, _label(pr.plan[pr.current][1]), message)
        flush(pr.io)
    end
    return pr
end

"""
    count_cycle!(pr)

Count a cycle simulated that is not reported on its own, such as a Newton
step tried and rejected, in the total at the end.
"""
count_cycle!(pr::CycleProgress) = (pr.cycles += 1; pr)

"""
    finish_phase!(pr, n, converged; note = "")

End the current phase after `n` cycles (or iterations): its line says whether
it converged, or that it is done if a later phase is planned.
"""
function finish_phase!(pr::CycleProgress, n, converged; note = "")
    k = pr.current
    value = isempty(pr.values) ? NaN : pr.values[end]
    next = k < length(pr.plan) ? pr.plan[k + 1][1] : nothing
    outcome = isnothing(next) ? (converged ? "converged" : "not converged") :
        (converged ? "done, on to $next" : "stopped, on to $next")
    parts = [_phase_columns(pr, k, "$(lpad(n, 3)) $(_plural(pr.noun, n))",
            "$(pr.measure) $(_sci(value))" * (isnan(pr.tol) ? "" : " $(value < pr.tol ? "<" : "≥") $(_sci(pr.tol))"),
            converged ? 1.0 : _log_fraction(pr.values, pr.tol)),
        outcome]
    isempty(note) || push!(parts, note)
    push!(parts, _duration(time() - pr.phase_started))
    _update_line!(pr, join(parts, "  "))
    isnothing(next) && _close_block!(pr)
    return pr
end

"""
    report_note!(pr, message)

Print a line that stays on screen, such as a heading or a change of plan,
below the phases so far.
"""
function report_note!(pr::CycleProgress, message)
    _close_block!(pr)
    println(pr.io, message)
    flush(pr.io)
    return pr
end

"""
    finish_progress!(pr, converged)

Print the closing line, with every cycle simulated and the total time.
"""
function finish_progress!(pr::CycleProgress, converged)
    word = converged ? "Converged" : "Not converged"
    report_note!(pr, "$word after $(pr.cycles) $(_plural("cycle", pr.cycles)) simulated in $(_duration(time() - pr.started))")
end

# ------------------------------------------------------------------------------
# Lines and the live block
# ------------------------------------------------------------------------------

const LABEL_WIDTH = 12
const COUNT_WIDTH = 14
const MEASURE_WIDTH = 25
const BAR_WIDTH = 12

_label(name) = "  " * rpad(name, LABEL_WIDTH)

# Where the text after the bar and its percentage starts
_status_column() = 2 + LABEL_WIDTH + 2 + COUNT_WIDTH + 2 + MEASURE_WIDTH + 2 + (BAR_WIDTH + 2) + 5 + 2

# The phase name, count, measure and bar, in columns so that the bars of the
# phases line up
function _phase_columns(pr, k, count, measure, frac)
    return _label(pr.plan[k][1]) * "  " * rpad(count, COUNT_WIDTH) * "  " * rpad(measure, MEASURE_WIDTH) * "  " *
        _bar(frac) * " " * lpad(round(Int, 100frac), 3) * "%"
end

function _waiting_line(pr, k)
    name, tol, noun = pr.plan[k]
    after = k > 1 ? "waiting for $(pr.plan[k - 1][1])" : "starting"
    return _phase_columns(pr, k, "", isnan(tol) ? "" : "to $(_sci(tol))", 0.0) * "  " * after
end

_plural(noun, n) = n == 1 ? noun : noun * "s"

# Set the line of the current phase: redraw the live block, or print it as a
# line of its own
function _update_line!(pr, line)
    pr.lines[pr.current] = line
    if pr.live
        _draw_block!(pr)
    else
        println(pr.io, line)
        flush(pr.io)
    end
    return pr
end

# Draw the lines of every planned phase over the ones drawn before, each kept
# within the terminal so that the block can be found again
function _draw_block!(pr)
    io = pr.io
    width = max(displaysize(io)[2] - 1, 20)
    buf = IOBuffer()
    print(buf, "\r")
    pr.drawn > 1 && print(buf, "\e[$(pr.drawn - 1)A")
    for (i, line) in enumerate(pr.lines)
        i > 1 && print(buf, "\n")
        print(buf, "\e[K", first(line, width))
    end
    write(io, take!(buf))
    flush(io)
    pr.drawn = length(pr.lines)
    pr.last_status = time()
    return nothing
end

# Leave the live block on screen, and print below it from now on
function _close_block!(pr)
    if pr.live && pr.drawn > 0
        println(pr.io)
        flush(pr.io)
    end
    pr.drawn = 0
    empty!(pr.plan)
    empty!(pr.lines)
    pr.current = 0
    return nothing
end

# ------------------------------------------------------------------------------
# Formatting
# ------------------------------------------------------------------------------

# Three significant figures in scientific notation, as 2.35e-03
function _sci(x)
    isfinite(x) || return string(x)
    x == 0 && return "0"
    e = floor(Int, log10(abs(x)))
    m = round(x / 10.0^e, digits = 2)
    abs(m) >= 10 && (m /= 10; e += 1)
    return string(m, "e", e < 0 ? "-" : "+", lpad(abs(e), 2, '0'))
end

function _duration(s)
    s < 10 && return "$(round(s, digits = 1)) s"
    s < 60 && return "$(round(Int, s)) s"
    s < 3600 && return "$(floor(Int, s / 60)) min $(round(Int, s % 60)) s"
    return "$(floor(Int, s / 3600)) h $(round(Int, (s % 3600) / 60)) min"
end

const BAR_EIGHTHS = collect("▏▎▍▌▋▊▉")

# A bar of `width` characters filled to `frac`, in eighths of a character
function _bar(frac; width = BAR_WIDTH)
    eighths = round(Int, clamp(frac, 0, 1) * 8width)
    full, part = divrem(eighths, 8)
    partial = part == 0 ? "" : string(BAR_EIGHTHS[part])
    return "▕" * "█"^full * partial * " "^(width - full - (part > 0)) * "▏"
end

# How far the smallest of `values` has come from the first towards `tol`, on
# a log scale
function _log_fraction(values, tol)
    first_, last_ = values[1], minimum(values)
    last_ < tol && return 1.0
    (first_ > tol && last_ > 0) || return 0.0
    return clamp(log(first_ / last_) / log(first_ / tol), 0, 1)
end

# Cycles still needed to bring `values` below `tol`, at the average rate of
# decrease of their smallest so far over the last `window` cycles, or nothing
# if that has not decreased. The first cycle, from the initial guess, is left
# out: its change is usually far larger than the rest.
function _cycles_to_reach(values, tol; window = 10)
    values = view(values, 2:length(values))
    n = length(values)
    n < 3 && return nothing
    m = min(window, n - 1)
    a, b = minimum(view(values, 1:n - m)), minimum(values)
    (a > 0 && b > 0) || return nothing
    b < tol && return 0
    rate = (b / a)^(1 / m)
    rate < 1 || return nothing
    return ceil(Int, log(tol / b) / log(rate))
end
