a = 0; b = 1
ϵ = 0.1
nₓ = 9; nₜ = 100; t_max = 1
left_bc = PeriodicBC(); right_bc = PeriodicBC()
u₀(x) = 1.

problem = AllenCahnProblem1D(a, b, ϵ, nₓ, nₜ, t_max, left_bc, right_bc, u₀)
assemble_system!(problem)
tsm = CrankNicolsonMethod()
u = solve(problem, tsm)

error = abs.(ones(length(u)) - u)
norm(error, Inf) < 1e-14