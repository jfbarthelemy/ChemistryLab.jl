using ChemistryLab
using Aqua
using DynamicQuantities
using ForwardDiff
using OptimaSolver
using OrdinaryDiffEq
using OrderedCollections
using Test
using TimerOutputs

macro testsection(str, block)
    return quote
        @timeit "$($(esc(str)))" begin
            @testset "$($(esc(str)))" begin
                $(esc(block))
            end
        end
    end
end

reset_timer!()

"""
    aqua_persistent_tasks(m::Module) -> Bool

`Aqua.test_persistent_tasks` in a child process with default bounds checking.

The check precompiles a synthetic package that depends on `m` and waits for a
sentinel written from that package's module body. Julia writes no
precompilation cache when `--check-bounds` is forced, so the body never runs,
the sentinel never appears, and the check fails for a reason that says nothing
about `m`. `Pkg.test()` forces the flag by default.

Running it in a child with `--check-bounds=auto` keeps both halves: the test
suite proper still runs with bounds checking on, and the check still runs.
"""
function aqua_persistent_tasks(m::Module)
    code = """
    using Aqua, $(nameof(m)), Test
    @testset "persistent_tasks" begin
        Aqua.test_persistent_tasks($(nameof(m)))
    end
    """
    cmd = `$(first(Base.julia_cmd())) --check-bounds=auto --startup-file=no --project=$(Base.active_project()) -e $code`
    return success(run(ignorestatus(cmd)))
end

@testsection "Aqua" begin
    Aqua.test_all(
        ChemistryLab;
        # `thermo_factories.jl` extends some thirty `Base` math functions to
        # `DynamicQuantities.Quantity`, stripping the unit before applying them.
        # That is type piracy by construction — the function and the type both
        # belong elsewhere — and it is deliberate: it is what lets a
        # dimensionless quantity be used wherever a number is expected. Declared
        # rather than hidden. It stays global, though: any code loading
        # ChemistryLab gets these methods whether it asked for them or not.
        piracies = (treat_as_own = [DynamicQuantities.Quantity],),
        # Run apart, see `aqua_persistent_tasks` above.
        persistent_tasks = false,
    )
    @test aqua_persistent_tasks(ChemistryLab)
end

@testsection "Construction tests" begin

    include("species.jl")
    include("formulas.jl")
    include("stoich_matrices.jl")
    include("databases.jl")
    include("parsing_utils.jl")
    include("reactions.jl")

end

@testsection "Speciation tests" begin
    include("speciation.jl")
end

@testsection "System and State tests" begin
    include("chemical_systems.jl")
    include("chemical_states.jl")
    include("test_volume_fractions.jl")
end

@testsection "Thermodynamics tests" begin
    include("thermodynamics.jl")
    include("hkf.jl")
end

@testsection "Equilibrium tests" begin
    include("activities.jl")
    include("solid_solutions.jl")
    include("equilibrium_reference.jl")
    include("published_values.jl")
    include("test_dual_solver.jl")
    include("certified_equilibrium.jl")
    include("equilibrium_constraints.jl")
end

@testsection "Utils tests" begin
    include("utils.jl")
end

@testsection "Kinetics tests" begin
    include("kinetics/test_rate_models.jl")
    include("kinetics/test_calorimetry.jl")
    include("kinetics/test_postprocessing.jl")
    include("kinetics/test_feasibility.jl")
    include("kinetics/test_calibration.jl")
    include("kinetics/test_implicit_step.jl")
    include("coupling_reference.jl")
end

print_timer()
println()
