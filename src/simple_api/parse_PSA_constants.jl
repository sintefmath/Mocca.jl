using JSON

function parse_PSA_constants(json_file::String)
    # Read and parse the JSON file
    json_data = JSON.parsefile(json_file)

    # Extract values from the JSON and initialize HaghpanahConstants
    constants = Mocca.PSAConstants{Float64}(
        molecularMassOfCO2 = json_data["properties"]["molecularMassOfCO2"]["value"],
        molecularMassOfN2 = json_data["properties"]["molecularMassOfN2"]["value"],
        R = json_data["properties"]["R"]["value"],
        Φ = json_data["properties"]["Φ"]["value"],
        b0 = SVector{2, Float64}(json_data["properties"]["b0"]["value"]...),
        d0 = SVector{2, Float64}(json_data["properties"]["d0"]["value"]...),
        ΔUbi = SVector{2, Float64}(json_data["properties"]["ΔUbi"]["value"]...),
        ΔUdi = SVector{2, Float64}(json_data["properties"]["ΔUdi"]["value"]...),
        qsbi = SVector{2, Float64}(json_data["properties"]["qsbi"]["value"]...),
        qsdi = SVector{2, Float64}(json_data["properties"]["qsdi"]["value"]...),
        ϵ_p = json_data["properties"]["ϵ_p"]["value"],
        D_m = json_data["properties"]["D_m"]["value"],
        τ = json_data["properties"]["τ"]["value"],
        d_p = json_data["properties"]["d_p"]["value"],
        V0_inter = json_data["properties"]["V0_inter"]["value"],
        fluid_viscosity = json_data["properties"]["fluid_viscosity"]["value"],
        K_z = json_data["properties"]["K_z"]["value"],
        K_w = json_data["properties"]["K_w"]["value"],
        ρ_s = json_data["properties"]["ρ_s"]["value"],
        ρ_g = json_data["properties"]["ρ_g"]["value"],
        C_pg = SVector{2, Float64}(json_data["properties"]["C_pg"]["value"]...),
        C_pa = SVector{2, Float64}(json_data["properties"]["C_pa"]["value"]...),
        C_ps = json_data["properties"]["C_ps"]["value"],
        r_in = json_data["properties"]["r_in"]["value"],
        r_out = json_data["properties"]["r_out"]["value"],
        h_in = json_data["properties"]["h_in"]["value"],
        h_out = json_data["properties"]["h_out"]["value"],
        ρ_w = json_data["properties"]["ρ_w"]["value"],
        C_pw = json_data["properties"]["C_pw"]["value"],
        T0 = json_data["properties"]["T0"]["value"],
        T_a = json_data["properties"]["T_a"]["value"],
        v_feed = json_data["properties"]["v_feed"]["value"],
        y_feed = SVector{2, Float64}(json_data["properties"]["y_feed"]["value"]...),
        p_high = json_data["properties"]["p_high"]["value"],
        p_intermediate = json_data["properties"]["p_intermediate"]["value"],
        p_low = json_data["properties"]["p_low"]["value"],
        λ = json_data["properties"]["λ"]["value"],
        T_feed = json_data["properties"]["T_feed"]["value"],
        L = json_data["properties"]["L"]["value"]
    )

    return constants
end