using JSON

function parse_PSA_constants(json_file::String)
    # Read and parse the JSON file
    json_data = JSON.parsefile(json_file)

    # Extract values from the JSON and initialize HaghpanahConstants
    constants = Mocca.PSAConstants{Float64}(
        # Physical constants
        molecularMassOfCO2 = json_data["physicalConstants"]["molecularMassOfCO2"]["value"],
        molecularMassOfN2 = json_data["physicalConstants"]["molecularMassOfN2"]["value"],
        R = json_data["physicalConstants"]["R"]["value"],
        # Dual-site Langmuir Isotherm
        b0 = SVector{2, Float64}(json_data["dslPars"]["b0"]["value"]...),
        d0 = SVector{2, Float64}(json_data["dslPars"]["d0"]["value"]...),
        ΔUbi = SVector{2, Float64}(json_data["dslPars"]["ΔUbi"]["value"]...),
        ΔUdi = SVector{2, Float64}(json_data["dslPars"]["ΔUdi"]["value"]...),
        qsbi = SVector{2, Float64}(json_data["dslPars"]["qsbi"]["value"]...),
        qsdi = SVector{2, Float64}(json_data["dslPars"]["qsdi"]["value"]...),
        # Adsorbent properties
        ϵ_p = json_data["adsorbentProps"]["ϵ_p"]["value"],
        D_m = json_data["adsorbentProps"]["D_m"]["value"],
        τ = json_data["adsorbentProps"]["τ"]["value"],
        d_p = json_data["adsorbentProps"]["d_p"]["value"],
        V0_inter = json_data["adsorbentProps"]["V0_inter"]["value"],
        ρ_s = json_data["adsorbentProps"]["ρ_s"]["value"],
        C_pa = SVector{2, Float64}(json_data["adsorbentProps"]["C_pa"]["value"]...),
        C_ps = json_data["adsorbentProps"]["C_ps"]["value"],
        # Column properties
        Φ = json_data["columnProps"]["Φ"]["value"],
        K_z = json_data["columnProps"]["K_z"]["value"],
        K_w = json_data["columnProps"]["K_w"]["value"],
        r_in = json_data["columnProps"]["r_in"]["value"],
        r_out = json_data["columnProps"]["r_out"]["value"],
        h_in = json_data["columnProps"]["h_in"]["value"],
        h_out = json_data["columnProps"]["h_out"]["value"],
        ρ_w = json_data["columnProps"]["ρ_w"]["value"],
        C_pw = json_data["columnProps"]["C_pw"]["value"],
        L = json_data["columnProps"]["L"]["value"],
        # Feed gas properties
        fluid_viscosity = json_data["feedProps"]["fluid_viscosity"]["value"],       
        ρ_g = json_data["feedProps"]["ρ_g"]["value"],
        C_pg = SVector{2, Float64}(json_data["feedProps"]["C_pg"]["value"]...),
        T_feed = json_data["feedProps"]["T_feed"]["value"],
        # Boundary conditions
        T_a = json_data["boundaryConditions"]["T_a"]["value"],
        v_feed = json_data["boundaryConditions"]["v_feed"]["value"],
        y_feed = SVector{2, Float64}(json_data["boundaryConditions"]["y_feed"]["value"]...),
        p_high = json_data["boundaryConditions"]["p_high"]["value"],
        p_intermediate = json_data["boundaryConditions"]["p_intermediate"]["value"],
        p_low = json_data["boundaryConditions"]["p_low"]["value"],
        λ = json_data["boundaryConditions"]["λ"]["value"],
        # Initial conditions
        T0 = json_data["initialConditions"]["T0"]["value"]
    )

    return constants
end