
abstract type AdsorptionSystem <: JutulDarcy.MultiComponentSystem end

const AdsorptionModel = Jutul.SimulationModel{<:Any,<:AdsorptionSystem,<:Any,<:Any}


# Overload JutulDarcy functions

JutulDarcy.number_of_components(sys::AdsorptionSystem) = sys.number_of_components
JutulDarcy.has_other_phase(::AdsorptionSystem) = false
JutulDarcy.phase_names(::AdsorptionSystem) = ["gas"]
JutulDarcy.number_of_phases(::AdsorptionSystem) = 1


# Specific functions needed for Mocca

function compute_permeability(p::ConstantsStruct)
    return 4 / 150 * ((p.Φ / (1 - p.Φ))^2) * (p.d_p / 2)^2 * p.Φ
end

function calc_dispersion(p::ConstantsStruct)
    return 0.7 * p.D_m + 0.5 * p.V0_inter * p.d_p
end


function compute_column_face_area(model::AdsorptionModel)\
    g = Jutul.physical_representation(model.data_domain)
    return g.deltas[2] * g.deltas[3]
end


function calc_bc_trans(model::AdsorptionModel, state)
    k = compute_permeability(model.system.p)
    dx = state.CellDx[1] / 2.0
    A = (π * model.system.p.r_in^2)
    return k * A / dx
end

function calc_bc_wall_trans(model::AdsorptionModel, state)
    k = model.system.p.K_w
    dx = state.CellDx[1] / 2.0
    A = area_wall(model.system)
    return k * A / dx
end


"Area of column wall [m^2]"
area_wall(sys::AdsorptionSystem) = (π * (sys.p.r_out^2 - sys.p.r_in^2))::Float64
area_wall_in(sys::AdsorptionSystem, Δx) = (π * sys.p.r_in * 2 * Δx)::Float64
area_wall_out(sys::AdsorptionSystem, Δx) = (π * sys.p.r_out * 2 * Δx)::Float64
