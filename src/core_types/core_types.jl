abstract type AdsorptionSystem <: JutulDarcy.MultiComponentSystem end

const AdsorptionModel = Jutul.SimulationModel{<:Any,<:AdsorptionSystem,<:Any,<:Any}

