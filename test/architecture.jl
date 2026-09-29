@testset "System hierarchy and aliases" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.FixedBed(constants)

    @test Mocca.AdsorptionSystem === Mocca.FixedBed
    @test Mocca.Column === Mocca.Unit
    @test system isa Mocca.SorbentBed
    @test system isa Mocca.DistributedUnit
    @test system isa Mocca.MoccaSystem

    model = Mocca.setup_process_model(system, constants; ncells = 5)
    @test model isa Mocca.FixedBedModel
    @test model isa Mocca.MoccaModel
end

@testset "Boundary quantities are parameters" begin
    constants = Mocca.HaghpanahConstants{Float64}()
    system = Mocca.FixedBed(constants)
    model = Mocca.setup_process_model(system, constants; ncells = 5)
    parameters = Mocca.setup_process_parameters(model)

    @test haskey(parameters, :Permeability)
    @test haskey(parameters, :BedCrossSectionArea)
    @test parameters[:BedCrossSectionArea][1] ≈ π * constants.r_in^2
    @test all(parameters[:Permeability] .≈ Mocca.compute_permeability(constants.Φ, constants.d_p))

    # The half-cell transmissibility follows the parameter, not the data domain
    state = merge(Dict{Symbol, Any}(:CellDx => parameters[:CellDx]), parameters)
    state[:Permeability] = 2 .* parameters[:Permeability]
    nt = NamedTuple(state)
    @test Mocca.calc_bc_trans(model, nt, 1) ≈ 2 * Mocca.calc_bc_trans(model, parameters |> NamedTuple, 1)
end
