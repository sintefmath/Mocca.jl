# ------------------------------------------------------------------------------
# Progress of the cyclic steady state solvers
# ------------------------------------------------------------------------------

"""
    CycleProgress(io = stdout; live = io isa Base.TTY)

Printout of the progress of the cyclic steady state solvers, given to them as
`progress`. Each solver is a phase (such as Anderson acceleration, then
Newton's method), and each cycle it simulates adds a line with the cycle count,
its measure of change against the tolerance, a bar of how far it has come (on
a log scale, from the first cycle to the tolerance, by the smallest change so
far, as the change of an accelerated cycle goes up and down), an estimate of
the cycles left from the recent rate of decrease of that smallest change, and
the time taken. A phase ends with a
line saying whether it converged.

With `live`, as in a terminal, the line of the current cycle is rewritten in
place, so only the end of each phase stays on screen.
"""
mutable struct CycleProgress
    io::IO
    live::Bool
    started::Float64
    # The current phase
    phase::String
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
    # When a live line was last written
    last_status::Float64
end

CycleProgress(io::IO = stdout; live = io isa Base.TTY) =
    CycleProgress(io, live, time(), "", "cycle", "change", NaN, true, time(), Float64[], Float64[], 0, 0.0)

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
start_phase!(::Nothing, args...; kwarg...) = nothing
report_cycle!(::Nothing, args...; kwarg...) = nothing
report_status!(::Nothing, args...; kwarg...) = nothing
finish_phase!(::Nothing, args...; kwarg...) = nothing
report_note!(::Nothing, args...) = nothing
finish_progress!(::Nothing, args...) = nothing

"""
    start_phase!(pr, name; tol, noun = "cycle", measure = "change", estimate = true)

Start a phase `name` that runs until its `measure` is below `tol`, counting
in `noun`s, with an estimate of the number left if `estimate`.
"""
function start_phase!(pr::CycleProgress, name; tol, noun = "cycle", measure = "change", estimate = true)
    pr.phase, pr.noun, pr.measure, pr.tol, pr.estimate = name, noun, measure, tol, estimate
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
    parts = [_phase_label(pr), "$(pr.noun) $(lpad(n, 3))",
        "$(pr.measure) $(_sci(value)) → $(_sci(pr.tol))", _bar(frac) * " " * lpad(round(Int, 100frac), 3) * "%"]
    isempty(note) || push!(parts, note)
    if !isnothing(left) && left > 0
        # The last few cycles, as the first ones include compilation
        recent = pr.durations[max(end - 4, 1):end]
        per_cycle = sum(recent) / length(recent)
        push!(parts, "~$left $(_plural(pr.noun, left)), ~$(_duration(left * per_cycle)) left")
    end
    push!(parts, _duration(time() - pr.phase_started))
    _print_line(pr, join(parts, "  "); persist = !pr.live)
    return pr
end

"""
    report_status!(pr, message; persist = false)

Show what the current phase is doing between cycles, such as building a
Jacobian. Unless `persist`, it is only shown on a live line, and at most four
times a second.
"""
function report_status!(pr::CycleProgress, message; persist = false)
    if persist
        _print_line(pr, _phase_label(pr) * "  " * message; persist = !pr.live)
    elseif pr.live && time() - pr.last_status > 0.25
        _print_line(pr, _phase_label(pr) * "  " * message; persist = false)
    end
    return pr
end

"""
    finish_phase!(pr, n, converged; note = "")

End the current phase after `n` cycles (or iterations), with a line that
stays on screen.
"""
function finish_phase!(pr::CycleProgress, n, converged; note = "")
    value = isempty(pr.values) ? NaN : pr.values[end]
    parts = [_phase_label(pr), "$(lpad(n, 3)) $(_plural(pr.noun, n))",
        "$(pr.measure) $(_sci(value)) $(value < pr.tol ? "<" : "≥") $(_sci(pr.tol))",
        converged ? "converged" : "not converged"]
    isempty(note) || push!(parts, note)
    push!(parts, _duration(time() - pr.phase_started))
    _print_line(pr, join(parts, "  "); persist = true)
    return pr
end

"""
    report_note!(pr, message)

Print a line that stays on screen, such as a heading or a change of plan.
"""
report_note!(pr::CycleProgress, message) = (_print_line(pr, message; persist = true); pr)

"""
    finish_progress!(pr, converged)

Print the closing line, with every cycle simulated and the total time.
"""
function finish_progress!(pr::CycleProgress, converged)
    word = converged ? "Converged" : "Not converged"
    report_note!(pr, "$word after $(pr.cycles) $(_plural("cycle", pr.cycles)) simulated in $(_duration(time() - pr.started))")
end

# ------------------------------------------------------------------------------
# Formatting
# ------------------------------------------------------------------------------

_phase_label(pr) = "  " * rpad(pr.phase, 12)

_plural(noun, n) = n == 1 ? noun : noun * "s"

function _print_line(pr, line; persist)
    io = pr.io
    if pr.live
        # Clear the live line, and keep a new one within the terminal so that
        # it can be rewritten
        print(io, "\r\e[K")
        persist || (line = first(line, max(displaysize(io)[2] - 1, 20)))
    end
    print(io, line)
    persist ? println(io) : (pr.last_status = time())
    flush(io)
    return nothing
end

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
function _bar(frac; width = 12)
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
