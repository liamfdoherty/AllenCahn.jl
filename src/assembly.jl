function assemble_system!(problem::AllenCahnProblem1D{TLBC, TRBC, TA}) where {TLBC <: NeumannBC, TRBC <: NeumannBC, TA <: SparseMatrixCSC}
    # TODO: Construct the system uₜ = -Au + f(u) to solve for u with timestepping
   
    # Construct first row corresponding to uₓ(a) = vₗ
    problem.A[1,1] = -1
    problem.A[1,2] = 1

    # Construct rows corresponding to interior mesh points
    for i ∈ 2:problem.nₓ - 1
        problem.A[i, i - 1] = 1
        problem.A[i,i] = -2
        problem.A[i, i + 1] = 1
    end

    # Construct last row corresponding to uₓ(b) = vᵣ (remember that the rhs has -vᵣ/Δx now since sign change for symmetry)
    problem.A[end, end] = -1
    problem.A[end, end - 1] = 1

    @. problem.A *= problem.ϵ/(problem.Δx)^2

    # TODO: Construct the right hand side of the system


    return problem
end