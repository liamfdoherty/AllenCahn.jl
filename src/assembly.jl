function assemble_system!(problem::AllenCahnProblem1D{TLBC, TRBC, TA}) where {TLBC <: NeumannBC, TRBC <: NeumannBC, TA <: SparseMatrixCSC}   
    # Construct first row corresponding to uₓ(a) = vₗ
    problem.A[1,1] = -2
    problem.A[1,2] = 2

    # Construct rows corresponding to interior mesh points
    for i ∈ 2:problem.nₓ + 1
        problem.A[i, i - 1] = 1
        problem.A[i,i] = -2
        problem.A[i, i + 1] = 1
    end

    # Construct last row corresponding to uₓ(b) = vᵣ
    problem.A[end, end] = -2
    problem.A[end, end - 1] = 2

    @. problem.A *= 1/(problem.Δx)^2
    
    return problem
end

function assemble_system!(problem::AllenCahnProblem1D{TLBC, TRBC, TA}) where {TLBC <: PeriodicBC, TRBC <: PeriodicBC, TA <: SparseMatrixCSC}   
    # Construct first row corresponding to u₀
    problem.A[1,1] = -2
    problem.A[1,2] = 1
    problem.A[1, end] = 1

    # Construct rows corresponding to interior mesh points
    for i ∈ 2:problem.nₓ
        problem.A[i, i - 1] = 1
        problem.A[i,i] = -2
        problem.A[i, i + 1] = 1
    end

    # Construct last row corresponding to uₙ
    problem.A[end, 1] = 1
    problem.A[end, end - 1] = 1
    problem.A[end, end] = -2

    @. problem.A *= 1/(problem.Δx)^2
    
    return problem
end