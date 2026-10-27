# Slope limiters for the advection of what moves through a unit with the
# sorbent or the heat transfer fluid (solid_advection.jl, heat_exchanger.jl).
#
# The value of an advected quantity φ at a face is reconstructed from the two
# cells upstream of it, as in MUSCL schemes:
#
#   φ_face = φ_up + ½·limited_slope(limiter, φ_up − φ_upup, φ_down − φ_up)
#
# The limiter keeps the reconstruction second order where the profile is
# smooth and falls back to first-order upwinding at extrema, so that no new
# peaks or dips appear (total variation diminishing, TVD). First-order
# upwinding smears fronts with a numerical diffusion of v·Δx/2.
#
# The limiters below are smooth in the slopes where these agree in sign. Where
# both slopes vanish they give 0, as first-order upwinding: their derivatives
# there are not defined, and a denominator that small would overflow them.
# The test is on a positive threshold, which compares values only; a test
# against 0 would be true for a dual number of value 0 with a positive
# derivative.

"""
    AbstractLimiter

How the value at a cell face of a quantity advected with the sorbent or the
heat transfer fluid is reconstructed from the cells upstream of it:

    φ_face = φ_up + ½·limited_slope(limiter, a, b)

where `a = φ_up − φ_upup` is the change across the upstream cell and
`b = φ_down − φ_up` the change downstream of it. [`Upwind`](@ref) is first
order; [`VanAlbada`](@ref), [`VanLeer`](@ref) and [`Minmod`](@ref) are second
order where the profile is smooth, and first order at extrema.

To add a limiter, define a subtype and a method of [`limited_slope`](@ref).
A limiter given as a function `ψ(r)` of the ratio of slopes `r = b/a` has
`limited_slope = ψ(b/a)·a`; avoid branches on the values where you can, as
the slope is differentiated in the Newton solver and in adjoints.
"""
abstract type AbstractLimiter end

"""
    limited_slope(limiter::AbstractLimiter, a, b) → s

The limited change `s` across the upstream cell, from the change across it,
`a`, and the change downstream, `b`; the face value is `φ_up + s/2`. A TVD
limiter gives `s = 0` where `a` and `b` differ in sign, and otherwise a value
between 0 and `2·min(|a|, |b|)` with the sign of `a`.
"""
function limited_slope end

"""
    Upwind()

First-order upwinding: the value at a face is that of the cell upstream of
it. Robust, but it smears fronts with a numerical diffusion of `v·Δx/2`.
"""
struct Upwind <: AbstractLimiter end
limited_slope(::Upwind, a, b) = zero(a * b)

"""
    VanAlbada()

The van Albada limiter, `ψ(r) = (r² + r)/(r² + 1)` for `r > 0` and 0
otherwise: `s = max(a·b, 0)·(a + b)/(a² + b²)`. Smooth where the profile is
monotone, and less compressive than [`VanLeer`](@ref).
"""
struct VanAlbada <: AbstractLimiter end
function limited_slope(::VanAlbada, a, b)
    ab = a * b
    d = a * a + b * b
    return d > _LIMITER_TINY^2 ? (ab + abs(ab)) / 2 * (a + b) / d : zero(ab)
end

"""
    VanLeer()

The van Leer limiter, `ψ(r) = (r + |r|)/(1 + |r|)`:
`s = 2·a·b/(a + b)` where `a` and `b` have the same sign, and 0 otherwise.
"""
struct VanLeer <: AbstractLimiter end
function limited_slope(::VanLeer, a, b)
    d = abs(a) + abs(b)
    return d > _LIMITER_TINY ? (a * abs(b) + abs(a) * b) / d : zero(a * b)
end

"""
    Minmod()

The minmod limiter, `ψ(r) = max(0, min(1, r))`: the smaller of `a` and `b` in
size where they have the same sign, and 0 otherwise. The most diffusive of
the TVD limiters.
"""
struct Minmod <: AbstractLimiter end
limited_slope(::Minmod, a, b) = (sign(a) + sign(b)) / 2 * min(abs(a), abs(b))

# Slopes below which a limiter gives 0, small enough not to matter and large
# enough that the square of a denominator does not underflow
const _LIMITER_TINY = 1e-75

"""
    advection_limiter(system) → AbstractLimiter

The limiter for what moves through a unit of `system` with its sorbent or its
heat transfer fluid. [`Upwind`](@ref) unless the system says otherwise.
"""
advection_limiter(system) = Upwind()

# Value at a face of a quantity advected from cell `up` towards cell `down`,
# with `upup` the cell upstream of `up`: φ(c) gives the value in cell c, and a
# missing neighbour (nothing) gives first-order upwinding there.
#
# Jutul finds the cells an equation depends on by tracing it once, and keeps
# only those (see networks/ports.jl). A limiter may not read a slope where it
# is flat, as at a uniform start, so both neighbours are also read with zero
# weight.
@inline function _face_value(limiter, φ, up, upup, down)
    φ_up = φ(up)
    (limiter isa Upwind || isnothing(upup) || isnothing(down)) && return φ_up
    φ_upup, φ_down = φ(upup), φ(down)
    return φ_up + limited_slope(limiter, φ_up - φ_upup, φ_down - φ_up) / 2 + _stencil_anchor(φ_upup) + _stencil_anchor(φ_down)
end

@inline _stencil_anchor(x) = zero(x) * x
