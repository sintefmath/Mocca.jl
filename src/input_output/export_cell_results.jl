using Mocca

function export_cell_results(outputfile::String, case::MoccaCase, states, timesteps; format="csv")

    if format != "csv"
        error("Currently only 'csv' format is supported.")
    end

    # Check if results file should be overwritten
    if isfile(outputfile)
        @error "Output file $outputfile already exists. Choose a new name or delete the existing file to save results."
        return
    end

    # Export to csv with columns: time, P, T, Tw, y1, ..., yn, q1, ..., qn
    open(outputfile, "w") do io
        # Write header
        header = ["time", "P", "T", "Tw"]
        for i in case.model.system.component_names
            push!(header, "y$(i)")
        end
        for i in case.model.system.component_names
            push!(header, "q$(i)")
        end
        println(io, join(header, ","))
        # Write data
        for (i, state) in enumerate(states)
            time = timesteps[i]
            P = state[:Pressure][end]
            T = state[:Temperature][end]
            Tw = state[:WallTemperature][end]
            y = state[:y][:,end]
            q = state[:AdsorbedConcentration][:,end]
            row = [time, P, T, Tw]
            append!(row, y)
            append!(row, q)
            println(io, join(row, ","))
        end
    end

    return
    
end
