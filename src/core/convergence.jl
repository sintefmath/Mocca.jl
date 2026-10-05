# Convergence criteria for the equations of distributed units, in measures that
# mean the same for any cell size and time step. They are checked against
# Jutul's `nonlinear_tolerance` (1e-3 by default):
#
# - CNV, for the gas and sorbent mass balances: the residual over a step as a
#   fraction of the cell's holdup, |R|·Δt/N. N is the gas in the cell for the
#   gas balance, and the gas and adsorbed moles for the loading equation.
# - ΔT [K], for the gas and wall energy balances: the residual over a step
#   divided by the heat capacity of the cell, the temperature error it amounts
#   to.
#
# The device equations have their own, in units/equipment/device.jl.

_cell_gas_moles(state, cell) = Jutul.value(state.TotalMolarConcentration[cell]) * Jutul.value(state.FluidVolume[cell])

function _cell_criterion(f, r, n_eq, dt)
    e = zeros(n_eq)
    for cell in axes(r, 2), i in 1:n_eq
        e[i] = max(e[i], abs(r[i, cell]) * dt / f(cell))
    end
    return e
end

_component_criterion_names(prefix, model) = map(n -> "$(prefix)_$n", component_names(model.system))

function Jutul.convergence_criterion(model::DistributedUnitModel, storage, eq::Jutul.ConservationLaw{:ComponentMasses}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    e = _cell_criterion(cell -> _cell_gas_moles(state, cell), r, size(r, 1), dt)
    return (CNV = (errors = e, names = _component_criterion_names("CNV", model)),)
end

function Jutul.convergence_criterion(model::SorbentBedModel, storage, eq::Jutul.ConservationLaw{:AdsorbedLoading}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    q = state.AdsorbedConcentration
    holdup(cell) = _cell_gas_moles(state, cell) + state.SolidVolume[cell] * sum(i -> Jutul.value(q[i, cell]), axes(q, 1))
    e = _cell_criterion(holdup, r, size(r, 1), dt)
    return (CNV = (errors = e, names = _component_criterion_names("CNV_q", model)),)
end

# Heat capacity [J/K] of a cell: the unit's heat capacity terms and the gas
function _cell_heat_capacity(model, state, cell)
    gas = _cell_gas_moles(state, cell) * Jutul.value(state.C_pg[cell]) * Jutul.value(state.AverageMolarMass[cell])
    return Jutul.value(sum_terms(heat_capacity, heat_capacity_terms(model.system), model, state, cell)) + gas
end

function Jutul.convergence_criterion(model::DistributedUnitModel, storage, eq::Jutul.ConservationLaw{:ColumnConservedEnergy}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    e = _cell_criterion(cell -> _cell_heat_capacity(model, state, cell), r, 1, dt)
    return (ΔT = (errors = e, names = ["ΔT"]),)
end

function Jutul.convergence_criterion(model::DistributedUnitModel, storage, eq::Jutul.ConservationLaw{:WallConservedEnergy}, eq_s, r;
        dt = 1.0, update_report = missing)
    state = storage.state
    wall(cell) = state.WallCrossSectionArea[1] * state.CellDx[cell] * state.WallDensity[1] * state.WallHeatCapacity[1]
    e = _cell_criterion(wall, r, 1, dt)
    return (ΔT = (errors = e, names = ["ΔT_wall"]),)
end
