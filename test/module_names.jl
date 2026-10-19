# Every global name that a method of Mocca uses must resolve in its module. A
# function body only looks its names up when it runs, so a name a module no
# longer sees, after code moves between Mocca's submodules, still loads; it
# fails later, when that line runs. This walks the lowered code of every method
# defined in Mocca's modules, including closures and keyword-argument bodies,
# and lists the names that do not resolve.

walk_globals!(acc, x::GlobalRef) = push!(acc, x)
walk_globals!(acc, x::Expr) = foreach(a -> walk_globals!(acc, a), x.args)
walk_globals!(acc, x) = nothing

function methods_defined_in(mods)
    ms = Method[]
    keep(m) = (m.module in mods && push!(ms, m); nothing)
    if isdefined(Core, :methodtable)
        # From Julia 1.12, every method is in one table
        Base.visit(keep, Core.methodtable)
    else
        # Before, each function and type has its own table
        tables = Set{Core.MethodTable}()
        for mod in union(Base.loaded_modules_array(), mods), n in names(mod; all = true, imported = true)
            isdefined(mod, n) || continue
            x = getfield(mod, n)
            T = Base.unwrap_unionall(x isa Type ? x : typeof(x))
            T isa DataType && isdefined(T.name, :mt) && push!(tables, T.name.mt)
        end
        foreach(mt -> Base.visit(keep, mt), tables)
    end
    return ms
end

# (module, name, where it is used) for each name used in `mods` that does not resolve
function unresolved_globals(mods)
    found = Set{Tuple{Module, Symbol, String}}()
    for m in methods_defined_in(mods)
        src = try
            Base.uncompressed_ast(m)
        catch
            continue
        end
        refs = GlobalRef[]
        foreach(st -> walk_globals!(refs, st), src.code)
        for g in refs
            if g.mod in mods && !isdefined(g.mod, g.name)
                push!(found, (g.mod, g.name, "$(m.file):$(m.line)"))
            end
        end
    end
    return found
end

function mocca_modules()
    mods = Set{Module}([Mocca])
    for n in names(Mocca; all = true)
        x = isdefined(Mocca, n) ? getfield(Mocca, n) : nothing
        x isa Module && x !== Mocca && parentmodule(x) === Mocca && push!(mods, x)
    end
    return mods
end

@testset "Every name used in Mocca's modules resolves" begin
    mods = mocca_modules()
    @test Mocca.Operation in mods && Mocca.Units in mods
    missing_names = unresolved_globals(mods)
    for (mod, name, location) in sort!(collect(missing_names); by = string)
        @error "$mod.$name is not defined, but is used at $location"
    end
    @test isempty(missing_names)
end

@testset "The check finds a name that does not resolve" begin
    planted = Module(:Planted)
    Core.eval(planted, :(in_body(x) = x + missing_in_body))
    Core.eval(planted, :(in_closure(x) = map(y -> y + missing_in_closure, x)))
    Core.eval(planted, :(in_keyword_body(x; k = 1) = x + k + missing_in_keyword_body))
    found = Set(r[2] for r in unresolved_globals(Set([planted])))
    @test found == Set([:missing_in_body, :missing_in_closure, :missing_in_keyword_body])
end
