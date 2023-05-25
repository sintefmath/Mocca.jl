using Parameters
import MAT

@with_kw struct AdsorptionParameters
    molecularMassOfCO2::Float64 = 44.01e-3 # kg / mole
    molecularMassOfN2::Float64 = 28e-3 # kg/mole
    R::Float64 = 8.3144598 # J⋅mol^−1⋅K^−1.
    Φ::Float64 = 0.37 # TODO: We should not hardcode this....
    b0::SVector{2,Float64} = @SVector [8.65e-7, 2.5e-6]
    d0::SVector{2,Float64} = @SVector [2.63e-8, 0.0]
    ΔUbi::SVector{2,Float64} = @SVector [-36_641.21, -1.58e4]
    ΔUdi::SVector{2,Float64} = @SVector [-35_690.66, 0.0]
    qsbi::SVector{2,Float64} = @SVector [3489.44, 6613.551]
    qsdi::SVector{2,Float64} = @SVector [2872.35, 0.00]
    ϵ_p::Float64 = 0.35
    D_m::Float64 = 1.6e-5
    τ::Float64 = 3.0

    "Particle diameter [m]"
    d_p::Float64 = 2e-3

    V0_inter::Float64 = 1.0

    fluid_viscosity::Float64 = 1.72e-5
    "[W/m/K]"
    K_z::Float64 = 0.0903 # TODO: Double check
    "[W/m/K]"
    K_w::Float64 = 16 # TODO: Triple check!!
    "Density of adsorbent, [kg m^{-3}]"
    ρ_s::Float64 = 1130
    "Specific heat capacity per component for fluid phase [J kg^{-1}K^{-1}]"
    C_pg::SVector{2,Float64} = @SVector [697.5687, 1096.4]
    "Specific heat capacity per component for adsorbent phase [J kg^{-1}K^{-1}]"
    C_pa::SVector{2,Float64} = @SVector [697.5687, 1096.4]
    "Specific heat capacity of solid adsorbent [J kg^[-1} K^{-1}]"
    C_ps::Float64 = 1070.0
    "Column radius [m]"
    r_in::Float64 = 1.445
    "Wall radius [m]"
    r_out::Float64 = 1.162
    "Heat transfer coefficient from column to wall [Wm^{-2}K^{-1}]" # TODO: Review this
    h_in::Float64 = 8.6
    "Heat transfer coefficient from wall to outside [Wm^{-2}K^{-1}]" # TODO: Review this
    h_out::Float64 = 2.5

    "Density of wall medium [kg m^{-3}]"
    ρ_w::Float64 = 7800.0 # TODO: Review this value and its documentation
    
    "Specific heat capacity for the wall [J kg^{-1}K^{-1}]"
    C_pw::Float64 = 502.0  # TODO: Review this value and its documentation

    "Initial temperature [K]"
    T0::Float64 = 298.15 # TODO: this value and its documentation
    "Ambient temperature [K]"
    T_a::Float64 = 298.15 # TODO: this value and its documentation

    # BC Stuff
    "Fluid velocity of feed gas [m s^{-1}]"
    v_feed::Float64 = 0.37

    "Mole fraction of the components [-]"
    y_feed::SVector{2, Float64} = [0.15, 0.85]
    
    "High pressure [Pa]"
    p_high::Float64 = 1e-5

    "Intermediate pressure [Pa]"
    p_intermediate::Float64 = 0.2e-5

    "Low pressure [Pa]"
    p_low::Float64 = 0.1e-5

    # TODO: Double check this
    "Pressure BC parameter [-]"
    λ::Float64 = 0.5

    "Feed temperature [K]"
    T_feed::Float64 = 298.15
end

function AdsorptionParameters(filename::String)
    data = MAT.matread(filename)

    if haskey(data, "problem")
        setup = data["problem"]["SimulatorSetup"]

    else
        setup = data["SimulatorSetup"]
    end
    model = setup["model"]
    schedule = setup["schedule"]
    bc = schedule["control"]["bc"]

    @show keys(model)
    rock = model["rock"]
    fluid = model["fluid"]
    separationsystem = model["separationSystem"]


#     face: 1
#     type: {'flux'}
#    value: @(t,PH,PI,PL,lambda)(PH-(PH-PL)*exp(-lambda*t))
#      sat: 1
#       Ta: 298.1500
#        y: [0.1500 0.8500]
#        T: 298.1500
# stepType: 'pressurisation'
#   stepNo: 1
#       PH: 100000
#       PI: 20000
#       PL: 10000
#   lambda: 0.5000
#    Vfeed: 0.3700
#  cycleNo: 1
#   stept0: 0

    parameters =  AdsorptionParameters(
        Φ = first(rock["poro"]),
        ρ_w = first(rock["rhoWall"]),
        ρ_s = first(rock["rhoAds"]),
        d_p = first(rock["rp"])*2,
        ϵ_p = first(rock["poroAds"]),
        τ = first(rock["tau"]),
        C_pw = first(rock["CpWall"]),
        C_ps = first(rock["CpAds"]),
        h_in = first(rock["hIn"]),
        h_out = first(rock["hOut"]),
        K_w = first(rock["lambdaWall"]),
        R = first(model["R"]),
        T_a = first(model["ambientTemperature"]),
        molecularMassOfCO2 = first(fluid["molarMass"]),
        molecularMassOfN2 = fluid["molarMass"][2],
        D_m = first(fluid["Dm"]),
        C_pa = SVector{2, Float64}(fluid["CpAds"]),
        C_pg = SVector{2, Float64}(fluid["CpAds"]),
        K_z = first(fluid["lambdaF"]),
        b0 = SVector{2, Float64}(separationsystem["b0"]),
        d0 = SVector{2, Float64}(separationsystem["d0"]),
        ΔUbi = SVector{2, Float64}(separationsystem["deltaUb"]),
        ΔUdi = SVector{2, Float64}(separationsystem["deltaUd"]),
        qsbi = SVector{2, Float64}(separationsystem["qsb"]),
        qsdi = SVector{2, Float64}(separationsystem["qsd"]),
        T0 = first(separationsystem["T0"]),
        y_feed = SVector{2, Float64}(bc["y"]),
        T_feed = first(bc["T"]),
        p_high = first(bc["PH"]),
        p_intermediate = first(bc["PI"]),
        p_low = first(bc["PL"]),
        λ = first(bc["lambda"]),
        v_feed = first(bc["Vfeed"])
        #ρ_g = first(rock[""]) # not used? TODO: Check
    )

    # Some sanity checks
    @assert parameters.D_m ≈ first(fluid["Dm"])
    @assert compute_permeability(parameters) ≈ first(rock["perm"])
    @assert axial_dispersion(parameters) ≈ first(fluid["DL"])

    return parameters
end

compute_permeability(p::AdsorptionParameters) = 4 / 150 * ((p.Φ / (1 - p.Φ))^2) * (p.d_p / 2)^2 * p.Φ

axial_dispersion(p::AdsorptionParameters) = 0.7 * p.D_m + 0.5 * p.V0_inter * p.d_p
