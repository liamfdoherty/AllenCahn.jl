function AllenCahnProblem1D(a, b, ϵ, nₓ, nₜ, left_bc::TLBC, right_bc::TRBC) where {TLBC <: NeumannBC, TRBC <: NeumannBC}
    @assert a < b "Cannot construct: Invalid interval!"
    @assert ϵ > 0 "Cannot construct: ϵ must be positive!"
    @assert (nₓ > 0 && nₜ > 0) "Cannot construct: nₓ and nₜ must both be positive!"

    Δt = (b - a)/nₜ
    return AllenCahnProblem1D(a, b, ϵ, nₓ, nₜ, Δt, left_bc, right_bc)
end