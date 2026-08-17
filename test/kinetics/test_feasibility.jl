using LinearAlgebra

# Helpers guarding the element balance of a re-speciation. They exist because a
# full OPC coupling stopped on a composition demanding 0.732 mol of sulfate
# against the 0.267 mol present — 174 % over — while the reported residual read
# 1.4e-2, because it was normalised by the 34 mol water budget.

@testset "_row_residual scales each element by its own budget" begin

    A = Float64[1 0; 0 1]

    # A violation on a small element must not be hidden by a large one. The
    # second element is out by 0.5, scaled by the larger of its budget (1.0) and
    # the matter flowing through the row (1.5), so a third — not 0.5/100, which
    # is what normalising by the largest total across all rows used to give.
    @test ChemistryLab._row_residual(A, [100.0, 1.5], [100.0, 1.0]) ≈ 1 / 3
    @test ChemistryLab._row_residual(A, [100.0, 1.5], [100.0, 1.0]) > 0.5 / 100

    # A row whose total is zero — electroneutrality — must not turn a rounding
    # error into a huge residual. Dividing by `b` alone reported 2.1e3 here.
    @test ChemistryLab._row_residual(A, [2.0, 1.0e-12], [2.0, 0.0]) < 1.0e-1

    # An exactly satisfied balance is exactly zero
    @test ChemistryLab._row_residual(A, [2.0, 3.0], [2.0, 3.0]) == 0.0

end

@testset "_budget_clip! caps a guess at the element budget" begin

    # One species consuming 3 units of an element that has 0.267 available: the
    # guess cannot exceed 0.089, whatever it was.
    n = [10.0]
    ChemistryLab._budget_clip!(n, Float64[3;;], [0.267])
    @test n[1] ≈ 0.267 / 3

    # A guess already within budget is untouched
    n2 = [0.01]
    ChemistryLab._budget_clip!(n2, Float64[3;;], [0.267])
    @test n2[1] ≈ 0.01

    # A negative total bounds nothing: H⁺ in a hydrating cement is met by the
    # hydroxides, and clipping on it would forbid every hydrate.
    n3 = [5.0]
    ChemistryLab._budget_clip!(n3, Float64[-2;;], [-11.8])
    @test n3[1] == 5.0

end

@testset "_restore_feasibility! lands in {A n = b, n ≥ 0}" begin

    A = Float64[1 1 0; 0 1 1]
    b = [1.0, 1.0]

    # Starting from a point that is both infeasible and negative
    n = [5.0, -3.0, 7.0]
    ChemistryLab._restore_feasibility!(n, A, b)
    @test maximum(abs, A * n .- b) < 1.0e-10
    @test all(>=(0), n)                        # ends on the box clamp, not the affine step

    # Already feasible: stays feasible, and cheaply
    n2 = [0.5, 0.5, 0.5]
    ChemistryLab._restore_feasibility!(n2, A, b)
    @test maximum(abs, A * n2 .- b) < 1.0e-10
    @test all(>=(0), n2)

    # The cement case that motivated it: ettringite carries 3 sulfate against a
    # 0.267 mol budget, so no feasible point may hold more than 0.089 of it.
    Ae = Float64[3 0; 2 2]                      # rows: sulfate, aluminium
    be = [0.267, 0.488]
    ne = [0.244, 0.0]                           # the point the coupling stalled on
    @test 3 * ne[1] > be[1]                     # infeasible as it stands
    ChemistryLab._restore_feasibility!(ne, Ae, be)
    @test maximum(abs, Ae * ne .- be) < 1.0e-10
    @test all(>=(0), ne)
    @test 3 * ne[1] <= be[1] + 1.0e-10          # sulfate no longer over-spent

end
