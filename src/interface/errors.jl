"""
    MoccaInputError(location, field, message)

An error in the description of a process, found when it is set up rather than
inside the solver. `location` says where it is, such as `"step 3 (blowdown) of
the cycle"` or `"adsorbent 13X"`, and `field` names the setting, or is empty.
[`validate`](@ref) returns them as a list instead of throwing.
"""
struct MoccaInputError <: Exception
    location::String
    field::String
    message::String
end
MoccaInputError(location, field::Symbol, message) = MoccaInputError(String(location), String(field), String(message))

function Base.showerror(io::IO, e::MoccaInputError)
    print(io, "MoccaInputError in ", e.location)
    isempty(e.field) || print(io, ", ", e.field)
    print(io, ": ", e.message)
end

input_error(location, field, message) = throw(MoccaInputError(String(location), String(field), String(message)))

# Run `f`, rethrowing a plain error as a MoccaInputError at `location`
function _at(f, location)
    try
        return f()
    catch e
        e isa MoccaInputError && rethrow()
        e isa InterruptException && rethrow()
        throw(MoccaInputError(String(location), "", sprint(showerror, e)))
    end
end
