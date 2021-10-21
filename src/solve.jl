function solve(problem::AllenCahnProblem1D{NeumannBC, NeumannBC, SparseArrays.SparseMatrixCSC{Float64, Int64}}, stepping_method::TTSM) where {TTSM <: BackwardEulerMethod}
    u = problem.u
    t = 0.
    t_trajectory = [t]; u_trajectory = [u]
    n = 1

    matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) - problem.Δt*problem.ϵ*problem.A
    boundary_correction = zeros(problem.nₓ + 2)
    boundary_correction[1] = -2*problem.ϵ*problem.left_bc.boundary_condition/problem.Δx
    boundary_correction[end] = 2*problem.ϵ*problem.right_bc.boundary_condition/problem.Δx

    while n <= problem.nₜ
        problem.rhs .-= problem.rhs
        problem.rhs .+= problem.f.(u)
        problem.rhs .+= boundary_correction
        b = u + problem.Δt*problem.rhs
        t += problem.Δt; push!(t_trajectory, t)
        u = matrix\b; push!(u_trajectory, u)
        n += 1
    end

    return t_trajectory, u_trajectory
end

function solve(problem::AllenCahnProblem1D{PeriodicBC, PeriodicBC, SparseArrays.SparseMatrixCSC{Float64, Int64}}, stepping_method::TTSM) where {TTSM <: BackwardEulerMethod}
    u = problem.u
    t = 0.
    t_trajectory = [t]; u_trajectory = [u]
    n = 1

    matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) - problem.Δt*problem.ϵ*problem.A

    while n <= problem.nₜ
        problem.rhs .-= problem.rhs
        problem.rhs .+= problem.f.(u)
        b = u + problem.Δt*problem.rhs
        t += problem.Δt; push!(t_trajectory, t)
        u = matrix\b; push!(u_trajectory, u)
        n += 1
    end

    push!(problem.x, problem.a + (problem.nₓ + 1)*problem.Δx)
    for i ∈ 1:length(u_trajectory)
        push!(u_trajectory[i], u_trajectory[i][1])
    end

    return t_trajectory, u_trajectory
end

function solve(problem::AllenCahnProblem1D{NeumannBC, NeumannBC, SparseArrays.SparseMatrixCSC{Float64, Int64}}, stepping_method::TTSM) where {TTSM <: CrankNicolsonMethod}
    u = problem.u
    t = 0.
    t_trajectory = [t]; u_trajectory = [u]
    n = 1

    left_matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) - (problem.Δt*problem.ϵ/2)*problem.A
    right_matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) + (problem.Δt*problem.ϵ/2)*problem.A
    
    while n <= problem.nₜ
        problem.rhs .-= problem.rhs
        problem.rhs .+= problem.f.(u)
        b = right_matrix*u + problem.Δt*problem.rhs
        t += problem.Δt; push!(t_trajectory, t)
        u = left_matrix\b; push!(u_trajectory, u)
        n += 1
    end
    
    return t_trajectory, u_trajectory
end

function solve(problem::AllenCahnProblem1D{PeriodicBC, PeriodicBC, SparseArrays.SparseMatrixCSC{Float64, Int64}}, stepping_method::TTSM) where {TTSM <: CrankNicolsonMethod}
    u = problem.u
    t = 0.
    t_trajectory = [t]; u_trajectory = [u]
    n = 1

    left_matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) - (problem.Δt*problem.ϵ/2)*problem.A
    right_matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) + (problem.Δt*problem.ϵ/2)*problem.A
    
    while n <= problem.nₜ
        problem.rhs .-= problem.rhs
        problem.rhs .+= problem.f.(u)
        b = right_matrix*u + problem.Δt*problem.rhs
        t += problem.Δt; push!(t_trajectory, t)
        u = left_matrix\b; push!(u_trajectory, u)
        n += 1
    end

    push!(problem.x, problem.a + (problem.nₓ + 1)*problem.Δx)
    for i ∈ 1:length(u_trajectory)
        push!(u_trajectory[i], u_trajectory[i][1])
    end

    return t_trajectory, u_trajectory
end