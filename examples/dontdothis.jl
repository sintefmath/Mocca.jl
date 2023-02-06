function update_our_total_masses!(
    totmass,
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    y,
    cTot,
    PhaseMassDensities,
    ix,
)
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:113 =#
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:123 =#
    sys = model.system
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:125 =#
    for cell in ix
        #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:127 =#
        totmass[1, cell] =
            cTot[cell] *
            y[1, cell] *
            PhaseMassDensities[1, cell] *
            model.domain.grid.pore_volumes[cell]
        #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:128 =#
        totmass[2, cell] =
            cTot[cell] *
            y[2, cell] *
            PhaseMassDensities[1, cell] *
            model.domain.grid.pore_volumes[cell]
        #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:131 =#
    end
end
#= /home/kjetil/projects/dac2023/Jutul.jl/src/variable_evaluation.jl:80 =#
function Jutul.get_dependencies(
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
)
    println("I am not going to be called")
    (:y, :cTot, :PhaseMassDensities)
end
#= /home/kjetil/projects/dac2023/Jutul.jl/src/variable_evaluation.jl:81 =#
function update_secondary_variable!(
    array_target,
    tv::JutulDarcy.TotalMasses,
    model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem},
    state,
    ix,
)
    update_our_total_masses!(
        array_target,
        tv,
        model,
        state.y,
        state.cTot,
        state.PhaseMassDensities,
        ix,
    )
end

#= /home/kjetil/projects/dac2023/Jutul.jl/src/variable_evaluation.jl:79 =#
function update_cTot!(
    ctot,
    tv::JutulDarcy.TotalMass,
    model::Jutul.SimulationModel{G,S},
    Pressure,
    Temperature,
    ix,
) where {G,S<:AdsorptionFlowSystem}
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:177 =#
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:179 =#
    sys = model.system
    #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:181 =#
    for cellindex in ix
        #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:182 =#
        ctot[cellindex] = Pressure[cellindex] / (sys.R * Temperature[cellindex])
        #= /home/kjetil/projects/dac2023/Mocca.jl/examples/macroexpansion.jl:185 =#
    end
end
#= /home/kjetil/projects/dac2023/Jutul.jl/src/variable_evaluation.jl:80 =#
function Jutul.get_dependencies(
    tv::JutulDarcy.TotalMass,
    model::Jutul.SimulationModel{G,S},
) where {G,S<:AdsorptionFlowSystem}
    println("I'm not going to be called either.")
    (:Pressure, :Temperature)
end
#= /home/kjetil/projects/dac2023/Jutul.jl/src/variable_evaluation.jl:81 =#
function update_secondary_variable!(
    array_target,
    tv::JutulDarcy.TotalMass,
    model::Jutul.SimulationModel{G,S},
    state,
    ix,
) where {G,S<:AdsorptionFlowSystem}
    update_cTot!(array_target, tv, model, state.Pressure, state.Temperature, ix)
end
