using Parameters
@with_kw struct AdsorptionFlowSystem <: JutulDarcy.MultiComponentSystem
    number_of_components::Int64 = 2
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
    d_p::Float64 = 2e-3
    V0_inter::Float64 = 0.03653 
    """
    Use this to turn on and off the forcing term in the equation.
    If you don't know what this does, leave it to 1.0.
    """
    forcing_term_coefficient::Float64 = 1.0
    fluid_viscosity::Float64 = 1.72e-5
    "[W/m/K]"
    K_z::Float64 = 0.0903 # TODO: Double check
    "Density of adsorbent, [kg m^{-3}]"
    ρ_s::Float64 = 1130
    "Specific heat capacity per component for fluid phase [J kg^{-1}K^{-1}]"
    C_pg::SVector{2, Float64} = @SVector [697.5687, 1096.4]
    "Specific heat capacity per component for adsorbent phase [J kg^{-1}K^{-1}]"
    C_pa::SVector{2, Float64} = @SVector [697.5687, 1096.4]
    "Specific heat capacity of solid adsorbent [J kg^[-1} K^{-1}]"
    C_ps::Float64 = 1070.
    "Column radius [m]"
    r_in::Float64 = 1.445
    "Wall radius [m]"
    r_out::Float64 = 1.162
    "Heat transfer coefficient from column to wall [Wm^{-2}K^{-1}]"
    h_in::Float64 = 8.6
end
const AdsorptionFlowModel = Jutul.SimulationModel{<:Any,<:AdsorptionFlowSystem,<:Any,<:Any}

JutulDarcy.number_of_components(sys::AdsorptionFlowSystem) = sys.number_of_components
JutulDarcy.has_other_phase(::AdsorptionFlowSystem) = false

compute_permeability(sys::AdsorptionFlowSystem) = 4 / 150 * ((sys.Φ / (1 - sys.Φ))^2) * (sys.d_p/2)^2 * sys.Φ
axial_dispersion(sys::AdsorptionFlowSystem) = 0.7 * sys.D_m + 0.5 * sys.V0_inter * sys.Φ * sys.d_p

"Area of column wall [m^2]"
area_wall(sys::AdsorptionFlowSystem) = π*(sys.r_out^2 - sys.r_in^2)