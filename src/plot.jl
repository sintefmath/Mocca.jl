import CairoMakie
import MakiePublication
import DelimitedFiles

# function get_pvar_symbols()
#     pvars = [:y, :Pressure, :AdsorbedConcentration, :Temperature, :WallTemperature]
#     return pvars
# end


# """
#     get_matlab_states(matfile)

# Convert matlab .mat file results struct to states like julia states struct
# """
# function get_matlab_states(inputfile)

#     pvars = get_pvar_symbols()
#     matfile = MAT.matread(inputfile)

#     results = matfile["results"]
#     key_to_file = Dict(
#         :Pressure => "pressure",
#         :Temperature => "T",
#         :WallTemperature => "Twall",
#         :y => "yCO2",
#         :AdsorbedConcentration => ["qCO2", "qN2"]
#     )
#     num_timesteps = size(results["pressure"])[2]

#     times = results["time"]

#     states = []
#     for t in 1:num_timesteps
#         tmp = Dict()
#         for (nsymb, symbol) in enumerate(pvars)
#             if symbol == :AdsorbedConcentration

#                 qCO2 = results[key_to_file[symbol][1]][:,t]
#                 qN2 = results[key_to_file[symbol][2]][:,t]
#                 matlabdata = hcat(qCO2,qN2)'
#             else
#                 k = key_to_file[symbol]
#                 matlabdata = results[k][:,t]

#             end
#             tmp[symbol] = matlabdata
#         end
#         push!(states,tmp)
#     end

#     return states, times

# end


# function get_vars_through_time(states,cell,varname)

#     return collect([states[i][varname][cell] for i in eachindex(states)])
# end

# function get_plot_from_symbol(x,y,varname)

# end




# function plot_pvars_outlet(model, all_sims)
#     return MakiePublication.with_theme(MakiePublication.theme_web()) do

#         pvars = get_pvar_symbols()
#         f = CairoMakie.Figure()

#         nc = size(all_sims[1][1][1][:Pressure], 1)
#         x = model.data_domain[:cell_centroids][1,:]

#         key_to_label = Dict(
#             :y => "y",
#             :Pressure => "p",
#             :AdsorbedConcentration => "q",
#             :Temperature => "T",
#             :WallTemperature => "T_{wall}"
#         )

#         outlet_cell = nc # TODO This should really be changed for EvacuationBC
#         for (nsymb, symbol) in enumerate(pvars)
#             @show symbol

#             if symbol == :y
#                 ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
#                 v = get_vars_through_time(all_sims[1][1],outlet_cell,symbol)
#                 times = Float32.(all_sims[1][2])
#                 CairoMakie.lines(ax, Float32.(times), Float32.(v))
#             end
#         end

#         CairoMakie.resize!(f.scene, (2 * 400, 3 * 400))
#         return f            

#     end
    
# end

# function plot_pvars_outlet(model, all_sims)
#     return MakiePublication.with_theme(MakiePublication.theme_web()) do

#         pvars = get_pvar_symbols()
#         f = CairoMakie.Figure()


        
        
#         nc = size(all_sims[1][1][1][:Pressure], 1)
#         x = model.data_domain[:cell_centroids][1,:]

#         key_to_label = Dict(
#             :y => "y",
#             :Pressure => "p",
#             :AdsorbedConcentration => "q",
#             :Temperature => "T",
#             :WallTemperature => "T_{wall}"
#         )

#         outlet_cell = nc # TODO This should really be changed for EvacuationBC
#         for (nsymb, symbol) in enumerate(pvars)
#             @show symbol
            
#             get_vars_through_time(states,cell,varname)

#             if symbol == :y
#                 ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
                
#                 for (states, times) in all_sims
#                     times = Float64.(times)
#                     v = collect([states[i][symbol][end] for i in eachindex(states)])
#                     CairoMakie.lines!(ax, times, v)
#                 end

#             elseif symbol == :AdsorbedConcentration
#                 for i in 1:size(all_sims[1][end][symbol], 1)
                    
#                     # Make axes
#                     ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), 
#                     xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")

#                     # Loop through simulation states
#                     for (states, times) in all_sims
#                         v = collect([states[j][symbol][end][i] for j in eachindex(states)])
#                         CairoMakie.lines!(ax, Float64.(times), Float64.(v))
#                     end
#                 end

#             else
#                 ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
#                 for (states, times) in all_sims
#                     v = collect([states[i][symbol][end] for i in eachindex(states)])
#                     CairoMakie.lines!(ax, Float64.(times), Float64.(v))
#                 end
#             end

#         end



#         CairoMakie.resize!(f.scene, (2 * 400, 3 * 400))
#         return f            

#     end
    
# end

function plot_outlet(model,states)
    return MakiePublication.with_theme(MakiePublication.theme_web()) do
        f = CairoMakie.Figure()
        nc = size(states[end][:Pressure], 1)
        x = model.data_domain[:cell_centroids][1,:]
        key_to_label = Dict(
            :y => "y",
            :Pressure => "p",
            :AdsorbedConcentration => "q",
            :Temperature => "T",
            :WallTemperature => "T_{wall}"
        )

        for (nsymb, symbol) in enumerate()
            @show symbol

            if size(states[end][symbol], 2) == 1
                ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
                CairoMakie.lines!(ax, t, Float64.([result[symbol][end] for result in states]), color=:darkgray)
            else
                for i in 1:size(states[end][symbol], 1)
                    ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")

                    CairoMakie.lines!(ax, t, Float64.([blah[symbol][i, end] for blah in states]), color=:darkgray)
                end
            end
        end
        CairoMakie.resize!(f.scene, (2 * 400, 3 * 400))
        return f
    end
end


function plot_states(model, states)
    return MakiePublication.with_theme(MakiePublication.theme_web()) do
        f = CairoMakie.Figure()
        nc = size(states[end][:Pressure], 1)
        x = model.data_domain[:cell_centroids][1,:]
        key_to_label = Dict(
            :y => "y",
            :Pressure => "p",
            :AdsorbedConcentration => "q",
            :Temperature => "T",
            :WallTemperature => "T_{wall}"
        )

        for (nsymb, symbol) in enumerate([:y, :Pressure, :AdsorbedConcentration, :Temperature, :WallTemperature])
            @show symbol
            if size(states[end][symbol], 2) == 1
                ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
                for j in eachindex(states)
                    CairoMakie.lines!(ax, x, Float64.(states[j][symbol][:]), color=:darkgray)
                end
                CairoMakie.lines!(ax, x, Float64.(states[end][symbol][:]), color=:red)

            else
                for i in 1:size(states[end][symbol], 1)
                    ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")
                    for j in eachindex(states)
                        CairoMakie.lines!(ax, x, Float64.(states[j][symbol][i, :]), color=:darkgray)
                    end

                    CairoMakie.lines!(ax, x, Float64.(states[end][symbol][i, :]), color=:red)

                end
            end
        end
        CairoMakie.resize!(f.scene, (2 * 400, 3 * 400))
        return f
    end
end






function plot_against_matlab_mat(states, inputfile, t::Float64, times_mocca::Vector{Float64})
    # Find corresponding index for Mocca

    index_mocca = argmin(abs.(times_mocca .- t))
    times_mrst = collect(Iterators.flatten(MAT.matread(inputfile)["results"]["time"]))
    index_mrst = argmin(abs.(times_mrst .- t))

    @info "Getting timesteps " index_mocca index_mrst times_mocca[index_mocca] times_mrst[index_mrst]

    return plot_against_matlab_mat(states, inputfile, timestep=index_mocca, timestep_matlab=index_mrst)
end

function plot_against_matlab_mat(states, inputfile; timestep=nothing, timestep_matlab=nothing)
    return MakiePublication.with_theme(MakiePublication.theme_web()) do

        f = CairoMakie.Figure()

        nc = size(states[end][:Pressure], 1)
        x = collect(LinRange(0.0, 1.0, nc))
        key_to_label = Dict(
            :y => "y",
            :Pressure => "p",
            :AdsorbedConcentration => "q",
            :Temperature => "T",
            :WallTemperature => "T_{wall}"
        )

        matfile = MAT.matread(inputfile)

        results = matfile["results"]
        key_to_file = Dict(
            :Pressure => "pressure",
            :Temperature => "T",
            :WallTemperature => "Twall",
            :y => ["yCO2"],
            :AdsorbedConcentration => ["qCO2", "qN2"]
        )
        if isnothing(timestep)
            timestep = size(states, 1)
        end

        timestep = min(timestep, size(states, 1))

        if isnothing(timestep_matlab)
            timestep_matlab = size(results[key_to_file[:Pressure]], 2)
        end

        for (nsymb, symbol) in enumerate([:y, :Pressure, :AdsorbedConcentration, :Temperature, :WallTemperature])
            @show symbol
            # Truncating to Float16 seems to be needed due to some weird cairomakie bug:
            # https://discourse.julialang.org/t/range-step-cannot-be-zero/66948/10
            # Reverting to Float64 so values can be compared with MRST. Not encountered bug so far.            
            # TODO: Fix the above
            if size(states[timestep][symbol], 2) == 1
                ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
                CairoMakie.lines!(ax, x, Float64.(states[timestep][symbol][:]), color=:red, label="Mocca.jl")

                if haskey(key_to_file, symbol)
                    k = key_to_file[symbol]

                    matlabdata = collect(Iterators.flatten(results[k][:, timestep_matlab]))

                    CairoMakie.lines!(ax, x, matlabdata, color=:grey, label="MRST")
                end
                CairoMakie.axislegend(ax)
            else
                for i in 1:size(states[timestep][symbol], 1)
                    ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")

                    CairoMakie.lines!(ax, x, Float64.(states[timestep][symbol][i, :]), color=:red, label="Mocca.jl")
                    if haskey(key_to_file, symbol)
                        keysforresult = key_to_file[symbol]

                        if size(keysforresult, 1) == 1
                            matlabdata = collect(Iterators.flatten(results[keysforresult[1]][:, timestep_matlab]))
                            if i > 1
                                matlabdata = 1.0 .- matlabdata
                            end
                        else
                            matlabdata = collect(Iterators.flatten(results[keysforresult[i]][:, timestep_matlab]))
                        end

                        CairoMakie.lines!(ax, x, matlabdata, color=:grey, label="MRST")
                        CairoMakie.axislegend(ax)
                    end
                end
            end

        end
        CairoMakie.resize!(f.scene, (2 * 400, 3 * 400))
        return f
    end

end

