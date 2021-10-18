function solve(problem::AllenCahnProblem1D, stepping_method::TTSM) where {TTSM <: BackwardEulerMethod}
    u = problem.u
    n = 1
    matrix = sparse(Matrix(1.0I, problem.nₓ + 2, problem.nₓ + 2)) - problem.Δt*problem.ϵ*problem.A
    while n <= problem.nₜ
        b = u + problem.Δt*problem.rhs
        u = matrix\b
        n += 1
    end
    return u
end