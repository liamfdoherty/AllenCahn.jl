"""
`AllenCahnProblem1D` - initialize the problem data structure
"""
function AllenCahnProblem1D(a, b, ϵ, nₓ, nₜ, t_max, left_bc::TLBC, right_bc::TRBC, u₀::TU) where {TLBC <: NeumannBC, TRBC <: NeumannBC, TU <: Function}
    @assert a < b "Cannot construct: Invalid interval!"
    @assert ϵ > 0 "Cannot construct: ϵ must be positive!"
    @assert (nₓ > 0 && nₜ > 0) "Cannot construct: nₓ and nₜ must both be positive!"

    Δx = (b - a)/(nₓ + 1)
    x = [a + i*Δx for i ∈ 0:nₓ + 1]
    Δt = t_max/nₜ
    f(x) = x - x^3
    u = u₀.(x)
    return AllenCahnProblem1D(a, b, ϵ, nₓ, Δx, x, nₜ, t_max, Δt, left_bc, right_bc, spzeros(nₓ + 2, nₓ + 2), zeros(nₓ + 2), f, u)
end