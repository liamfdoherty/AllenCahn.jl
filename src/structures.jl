"""
`AllenCahnProblem1D` - structure to solve the one-dimensional Allen-Cahn equation, that is, solves

uₜ = -Auₓₓ + f(u)

on (a, b), where A = -Δ and f(u) = -u³ + u.
"""
struct AllenCahnProblem1D{TLBC, TRBC, TA}
    a::Real # Left endpoint
    b::Real # Right endpoint
    ϵ::Real # Model parameter
    nₓ::Int # Total number of spatial mesh points
    Δx::Real # Spatial mesh spacing
    x::Vector # Interior of the spatial mesh
    nₜ::Int # Number of timesteps
    Δt::Real # Length of timestep
    left_bc::TLBC # Boundary condition on the left endpoint
    right_bc::TRBC # Boundary condition on the right endpoint
    A::TA # Matrix approximation of differential operator for the problem
    rhs::Vector # Right hand side of the end system
end

"""
`NeumannBC` - structure to house a Neumann boundary condition, i.e., a boundary condition for which

uₓ(a) = vₗ  or  uₓ(b) = vᵣ
"""
struct NeumannBC
    boundary_condition::Real
end