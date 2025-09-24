using JSON

function parse_PSA_constants(json_file::String)
    # Read and parse the JSON file
    json_data = JSON.parsefile(json_file)

    # Extract values from the JSON and initialize HaghpanahConstants
    constants = Mocca.PSAConstants{Float64}(
        molecularMassOfCO2 = json_data["molecularMassOfCO2"]["value"],
        molecularMassOfN2 = json_data["molecularMassOfN2"]["value"],
        R = json_data["R"]["value"],
        Φ = json_data["Φ"]["value"],
        b0 = SVector{2, Float64}(json_data["b0"]["value"]...),
        d0 = SVector{2, Float64}(json_data["d0"]["value"]...),
        ΔUbi = SVector{2, Float64}(json_data["ΔUbi"]["value"]...),
        ΔUdi = SVector{2, Float64}(json_data["ΔUdi"]["value"]...),
        qsbi = SVector{2, Float64}(json_data["qsbi"]["value"]...),
        qsdi = SVector{2, Float64}(json_data["qsdi"]["value"]...),
        ϵ_p = json_data["ϵ_p"]["value"],
        D_m = json_data["D_m"]["value"],
        τ = json_data["τ"]["value"],
        d_p = json_data["d_p"]["value"],
        V0_inter = json_data["V0_inter"]["value"],
        fluid_viscosity = json_data["fluid_viscosity"]["value"],
        K_z = json_data["K_z"]["value"],
        K_w = json_data["K_w"]["value"],
        ρ_s = json_data["ρ_s"]["value"],
        ρ_g = json_data["ρ_g"]["value"],
        C_pg = SVector{2, Float64}(json_data["C_pg"]["value"]...),
        C_pa = SVector{2, Float64}(json_data["C_pa"]["value"]...),
        C_ps = json_data["C_ps"]["value"],
        r_in = json_data["r_in"]["value"],
        r_out = json_data["r_out"]["value"],
        h_in = json_data["h_in"]["value"],
        h_out = json_data["h_out"]["value"],
        ρ_w = json_data["ρ_w"]["value"],
        C_pw = json_data["C_pw"]["value"],
        T0 = json_data["T0"]["value"],
        T_a = json_data["T_a"]["value"],
        v_feed = json_data["v_feed"]["value"],
        y_feed = SVector{2, Float64}(json_data["y_feed"]["value"]...),
        p_high = json_data["p_high"]["value"],
        p_intermediate = json_data["p_intermediate"]["value"],
        p_low = json_data["p_low"]["value"],
        λ = json_data["λ"]["value"],
        T_feed = json_data["T_feed"]["value"],
        L = json_data["L"]["value"]
    )

    return constants
end