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

    # Construct the right hand side of the system (nonlinear term plus correction for nonhomogeneous boundary conditions)
    problem.rhs .+= problem.f.(problem.u)
    boundary_correction = zeros(problem.nₓ + 2)
    boundary_correction[1] = -2*problem.left_bc.boundary_condition/problem.Δx
    boundary_correction[end] = 2*problem.right_bc.boundary_condition/problem.Δx
    problem.rhs .+= boundary_correction
    
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

    # Construct the right hand side of the system (nonlinear term, there is no correction for BCs)
    problem.rhs .+= problem.f.(problem.u)
    
    return problem
end