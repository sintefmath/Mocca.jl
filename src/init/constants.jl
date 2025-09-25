using Parameters
abstract type ConstantsStruct end

@with_kw struct HaghpanahConstants{RealT} <: ConstantsStruct
    molecularMassOfCO2::RealT = 44.01e-3 # kg / mole
    molecularMassOfN2::RealT = 28e-3 # kg/mole
    R::RealT = 8.3144598 # J⋅mol^−1⋅K^−1. 
    Φ::RealT = 0.37 # TODO: We should not hardcode this....
    b0::SVector{2,RealT} = @SVector [8.65e-7, 2.5e-6]
    d0::SVector{2,RealT} = @SVector [2.63e-8, 0.0]
    ΔUbi::SVector{2,RealT} = @SVector [-36_641.21, -1.58e4]
    ΔUdi::SVector{2,RealT} = @SVector [-35_690.66, 0.0]
    qsbi::SVector{2,RealT} = @SVector [3489.44, 6613.551]
    qsdi::SVector{2,RealT} = @SVector [2872.35, 0.00]
    ϵ_p::RealT = 0.35
    D_m::RealT = 1.6e-5
    τ::RealT = 3.0

    "Particle diameter [m]"
    d_p::RealT = 2e-3

    V0_inter::RealT = 1.0

    fluid_viscosity::RealT = 1.72e-5
    "[W/m/K]"
    K_z::RealT = 0.0903 # Haghpanah et al. use 0.09
    "[W/m/K]"
    K_w::RealT = 16.0
    "Density of adsorbent, [kg m^{-3}]"
    ρ_s::RealT = 1130.0
    "Density of gas, [kg m^{-3}]"
    ρ_g::RealT = 1.22638310956

    "Specific heat capacity per component for fluid phase [J kg^{-1}K^{-1}]"
    C_pg::SVector{2,RealT} = @SVector [697.5687, 1096.4]
    "Specific heat capacity per component for adsorbent phase [J kg^{-1}K^{-1}]"
    C_pa::SVector{2,RealT} = @SVector [697.5687, 1096.4]
    "Specific heat capacity of solid adsorbent [J kg^[-1} K^{-1}]"
    C_ps::RealT = 1070.0
    "Column radius [m]"
    r_in::RealT = 0.1445
    "Wall radius [m]"
    r_out::RealT = 0.162
    "Heat transfer coefficient from column to wall [Wm^{-2}K^{-1}]"
    h_in::RealT = 8.6
    "Heat transfer coefficient from wall to outside [Wm^{-2}K^{-1}]"
    h_out::RealT = 2.5

    "Density of wall medium [kg m^{-3}]"
    ρ_w::RealT = 7800.0

    "Specific heat capacity for the wall [J kg^{-1}K^{-1}]"
    C_pw::RealT = 502.0

    "Initial temperature [K]"
    T0::RealT = 298.15
    "Ambient temperature [K]"
    T_a::RealT = 298.15

    # BC Stuff
    "Fluid velocity of feed gas [m s^{-1}]"
    v_feed::RealT = 0.37

    "Mole fraction of the components [-]"
    y_feed::SVector{2,RealT} = [0.15, 0.85]

    "p_high High pressure [Pa]"
    p_high::RealT = 1e5

    "Intermediate pressure [Pa]"
    p_intermediate::RealT = 0.2e5

    "Low pressure [Pa]"
    p_low::RealT = 0.1e5

    "Pressure BC parameter [s^{-1}]"
    λ::RealT = 0.5

    "Feed temperature [K]"
    T_feed::RealT = 298.15

    "Column length [m]"
    L::RealT = 1.0

end

@with_kw mutable struct PSAConstants{RealT} <: ConstantsStruct
    molecularMassOfCO2::RealT
    molecularMassOfN2::RealT
    R::RealT
    Φ::RealT
    b0::SVector{2,RealT}
    d0::SVector{2,RealT}
    ΔUbi::SVector{2,RealT}
    ΔUdi::SVector{2,RealT}
    qsbi::SVector{2,RealT}
    qsdi::SVector{2,RealT}
    ϵ_p::RealT
    D_m::RealT
    τ::RealT

    "Particle diameter [m]"
    d_p::RealT

    V0_inter::RealT

    fluid_viscosity::RealT
    "[W/m/K]"
    K_z::RealT
    "[W/m/K]"
    K_w::RealT
    "Density of adsorbent, [kg m^{-3}]"
    ρ_s::RealT
    "Density of gas, [kg m^{-3}]"
    ρ_g::RealT

    "Specific heat capacity per component for fluid phase [J kg^{-1}K^{-1}]"
    C_pg::SVector{2,RealT}
    "Specific heat capacity per component for adsorbent phase [J kg^{-1}K^{-1}]"
    C_pa::SVector{2,RealT}
    "Specific heat capacity of solid adsorbent [J kg^[-1} K^{-1}]"
    C_ps::RealT
    "Column radius [m]"
    r_in::RealT
    "Wall radius [m]"
    r_out::RealT
    "Heat transfer coefficient from column to wall [Wm^{-2}K^{-1}]"
    h_in::RealT
    "Heat transfer coefficient from wall to outside [Wm^{-2}K^{-1}]"
    h_out::RealT

    "Density of wall medium [kg m^{-3}]"
    ρ_w::RealT

    "Specific heat capacity for the wall [J kg^{-1}K^{-1}]"
    C_pw::RealT

    "Initial temperature [K]"
    T0::RealT
    "Ambient temperature [K]"
    T_a::RealT

    # BC Stuff
    "Fluid velocity of feed gas [m s^{-1}]"
    v_feed::RealT

    "Mole fraction of the components [-]"
    y_feed::SVector{2,RealT}

    "p_high High pressure [Pa]"
    p_high::RealT

    "Intermediate pressure [Pa]"
    p_intermediate::RealT

    "Low pressure [Pa]"
    p_low::RealT

    "Pressure BC parameter [s^{-1}]"
    λ::RealT

    "Feed temperature [K]"
    T_feed::RealT

    "Column length [m]"
    L::RealT

end


#function HaghpanahConstants(; kwarg...)
#    HaghpanahConstants{Float64}(; kwarg...)
#end


# axial_dispersion(p::HaghpanahConstants) = 0.7 * p.D_m + 0.5 * p.V0_inter * p.d_p
