using Mocca
using JSON3

function parse_input(input_dict::Dict{String, Any})
        return Mocca.setup_case_from_dict(input_dict)
end

function parse_input(filepath::String)

    input_dict = JSON3.read(open(filepath), Dict{String, Any})
    adsorptionInput = parse_input_from_dict(input_dict)
    
    return adsorptionInput
end


function parse_input_from_dict(input_dict::Dict{String, Any})

   # Extract values from the JSON and initialize HaghpanahConstants
    constants = Mocca.adsorptionInput{Float64}(
        # Physical constants
        molecularMassOfCO2 = input_dict["physicalConstants"]["molecularMassOfCO2"]["value"],
        molecularMassOfN2 = input_dict["physicalConstants"]["molecularMassOfN2"]["value"],
        R = input_dict["physicalConstants"]["R"]["value"],
        # Dual-site Langmuir Isotherm
        b0 = SVector{2, Float64}(input_dict["dslPars"]["b0"]["value"]...),
        d0 = SVector{2, Float64}(input_dict["dslPars"]["d0"]["value"]...),
        ΔUbi = SVector{2, Float64}(input_dict["dslPars"]["ΔUbi"]["value"]...),
        ΔUdi = SVector{2, Float64}(input_dict["dslPars"]["ΔUdi"]["value"]...),
        qsbi = SVector{2, Float64}(input_dict["dslPars"]["qsbi"]["value"]...),
        qsdi = SVector{2, Float64}(input_dict["dslPars"]["qsdi"]["value"]...),
        # Adsorbent properties
        ϵ_p = input_dict["adsorbentProps"]["ϵ_p"]["value"],
        D_m = input_dict["adsorbentProps"]["D_m"]["value"],
        τ = input_dict["adsorbentProps"]["τ"]["value"],
        d_p = input_dict["adsorbentProps"]["d_p"]["value"],
        V0_inter = input_dict["adsorbentProps"]["V0_inter"]["value"],
        ρ_s = input_dict["adsorbentProps"]["ρ_s"]["value"],
        C_pa = SVector{2, Float64}(input_dict["adsorbentProps"]["C_pa"]["value"]...),
        C_ps = input_dict["adsorbentProps"]["C_ps"]["value"],
        # Column properties
        Φ = input_dict["columnProps"]["Φ"]["value"],
        K_z = input_dict["columnProps"]["K_z"]["value"],
        K_w = input_dict["columnProps"]["K_w"]["value"],
        r_in = input_dict["columnProps"]["r_in"]["value"],
        r_out = input_dict["columnProps"]["r_out"]["value"],
        h_in = input_dict["columnProps"]["h_in"]["value"],
        h_out = input_dict["columnProps"]["h_out"]["value"],
        ρ_w = input_dict["columnProps"]["ρ_w"]["value"],
        C_pw = input_dict["columnProps"]["C_pw"]["value"],
        L = input_dict["columnProps"]["L"]["value"],
        # Feed gas properties
        fluid_viscosity = input_dict["feedProps"]["fluid_viscosity"]["value"],       
        ρ_g = input_dict["feedProps"]["ρ_g"]["value"],
        C_pg = SVector{2, Float64}(input_dict["feedProps"]["C_pg"]["value"]...),
        T_feed = input_dict["feedProps"]["T_feed"]["value"],
        # Boundary conditions
        T_a = input_dict["boundaryConditions"]["T_a"]["value"],
        v_feed = input_dict["boundaryConditions"]["v_feed"]["value"],
        y_feed = SVector{2, Float64}(input_dict["boundaryConditions"]["y_feed"]["value"]...),
        p_high = input_dict["boundaryConditions"]["p_high"]["value"],
        p_intermediate = input_dict["boundaryConditions"]["p_intermediate"]["value"],
        p_low = input_dict["boundaryConditions"]["p_low"]["value"],
        λ = input_dict["boundaryConditions"]["λ"]["value"],
        # Initial conditions
        P_init = input_dict["initialConditions"]["P_init"]["value"],
        T0 = input_dict["initialConditions"]["T0"]["value"],
        Tw_init = input_dict["initialConditions"]["Tw_init"]["value"],
        y_init = SVector{2, Float64}(input_dict["initialConditions"]["y_init"]["value"]...),

        # Process specification
        stage_types = input_dict["processSpecification"]["stage_types"]["value"],
        stage_durations = input_dict["processSpecification"]["stage_durations"]["value"],
        num_cycles = input_dict["processSpecification"]["num_cycles"]["value"],

        # Simulation parameters
        system_type = input_dict["simulation"]["system_type"]["value"],
        ncells = input_dict["simulation"]["ncells"]["value"],
        maxdt = input_dict["simulation"]["maxdt"]["value"],
        timestep_selectors = input_dict["simulation"]["timestep_selectors"]["value"],
        # Solver parameters
        linear_solver = input_dict["solver"]["linear_solver"]["value"],
        info_level = input_dict["solver"]["info_level"]["value"]
    )

    return constants
end