# Gradient checks of every process of the library, and of the processes the
# examples build, against central differences (see check_gradients). They take
# about a quarter of an hour, so they run with the tests only when
# MOCCA_GRADIENT_CHECKS is set, and MOCCA_GRADIENT_CHECKS_ONLY picks cases by
# part of their names:
#
#     MOCCA_GRADIENT_CHECKS=true julia --project -e 'using Pkg; Pkg.test()'
#     MOCCA_GRADIENT_CHECKS=true MOCCA_GRADIENT_CHECKS_ONLY=MBTSA,LPP julia ...
#
# The examples build these processes: the four-stage VSA (getting started,
# process from parts, gradients, optimisation, flowsheet by hand), the two-bed
# VSA, the wet flue gas LPP and dual-adsorbent VSAs and the MBTSA (moving bed
# and its optimisation). The breakthrough and history matching examples run a
# single step once, not to a steady state; history matching takes its
# gradients from Jutul's adjoints directly.
using Test
using Mocca
using Mocca.Units

function check(p)
    r = check_gradients(p)
    @test r.passed
    @test count(row -> row.status == :pass, r.rows) >= length(r.rows) / 2
    return r
end

# Only the cases named in MOCCA_GRADIENT_CHECKS_ONLY, comma-separated, if set
_only = split(get(ENV, "MOCCA_GRADIENT_CHECKS_ONLY", ""), ",", keepempty = false)
run_case(name) = isempty(_only) || any(k -> occursin(k, name), _only)

@testset "Gradient checks" begin
    run_case("four-stage VSA") && @testset "four-stage VSA" begin
        # The file's own variables: the feed velocity, the blowdown and
        # evacuation pressures and the durations of three steps
        p = read_process("haghpanah_2013_four_stage_vsa")
        add_variable!(p, :length, (:geometry, :Bed, :length))
        add_variable!(p, :density, (:parameter, :Bed, :AdsorbentDensity))
        check(p)
    end
    run_case("two-bed VSA with equalisation") && @testset "two-bed VSA with equalisation" begin
        p = two_bed_vsa_process()
        add_variable!(p, :v_feed, ("A adsorption", :velocity), ("B adsorption", :velocity))
        add_variable!(p, :P_L, ((k, :to) for k in (2, 4, 7, 9))...)
        add_variable!(p, :t_eq, (5, :duration), (10, :duration))
        add_variable!(p, :C_eq, (5, :conductance), (10, :conductance))
        check(p)
    end
    run_case("LPP VSA, wet flue gas") && @testset "LPP VSA, wet flue gas" begin
        p = read_process("krishnamurthy_2014_lpp_vsa")
        add_variable!(p, :v_feed, ("adsorption", :velocity))
        add_variable!(p, :P_L, ("evacuation", :to))
        add_variable!(p, :t_ads, ("adsorption", :duration))
        add_variable!(p, :length, (:geometry, only(Mocca.beds(p.flowsheet)), :length))
        check(p)
    end
    run_case("dual-adsorbent VSA, wet flue gas") && @testset "dual-adsorbent VSA, wet flue gas" begin
        p = read_process("krishnamurthy_2014_dual_adsorbent_vsa")
        add_variable!(p, :v_feed, ("silica gel", "adsorption", :velocity))
        add_variable!(p, :P_L_13X, ("13X", "evacuation", :to))
        add_variable!(p, :t_ads, ("silica gel", "adsorption", :duration), ("13X", "adsorption", :duration))
        add_variable!(p, :silica_gel_length, (:geometry, :SilicaGel, :length))
        check(p)
    end
    run_case("MBTSA") && @testset "MBTSA" begin
        p = read_process("mondino_2022_mbtsa")
        sections = (:Adsorption, :Preheating, :Desorption, :Precooling, :Cooling)
        add_variable!(p, :purge, (:setting, :V_purge, :rate))
        add_variable!(p, :hot_water, (:parameter, :Desorption, :HXFluidInletTemperature))
        add_variable!(p, :sorbent_flow, ((:parameter, s, :SorbentMassFlow) for s in sections)...)
        add_variable!(p, :desorption_height, (:geometry, :Desorption, :length))
        add_variable!(p, :extraction, (:setting, :V_product, :outlet, :pressure))
        check(p)
    end
end
