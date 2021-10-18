function solve(problem::AllenCahnProblem1D, stepping_method::TTSM) where {TTSM <: BackwardEulerMethod}
    u = problem.u
    n = 1
    matrix = sparse(Matrix(1.0I, size(problem.A)[1], size(problem.A)[2])) - problem.Δt*problem.ϵ*problem.A
    while n <= problem.nₜ
        b = u + problem.Δt*problem.rhs
        u = matrix\b
        n += 1
    end
    if typeof(problem.left_bc) <: PeriodicBC && typeof(problem.right_bc) <: PeriodicBC
        # Add the last point
        push!(problem.x, problem.a + (problem.nₓ + 1)*problem.Δx)
        push!(u, u[1])
    end
    return u
end