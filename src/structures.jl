"""
`AllenCahnProblem1D` - structure to solve the one-dimensional Allen-Cahn equation, that is, solves

uₜ = Auₓₓ + f(u)

on (a, b), where A = Δ and f(u) = -u³ + u.
"""
struct AllenCahnProblem1D{TLBC, TRBC, TA}
    a::Real # Left endpoint
    b::Real # Right endpoint
    ϵ::Real # Model parameter
    nₓ::Int # Total number of interior spatial mesh points
    Δx::Real # Spatial mesh spacing
    x::Vector # Spatial mesh points including boundary
    nₜ::Int # Number of timesteps
    t_max::Real # Final time for integration
    Δt::Real # Length of timestep
    left_bc::TLBC # Boundary condition on the left endpoint
    right_bc::TRBC # Boundary condition on the right endpoint
    A::TA # Matrix approximation of differential operator for the problem
    rhs::Vector # Right hand side of the end system
    f::Function # Nonlinear part of semilinear ODE system (u - u³ for this problem)
    u::Vector # Solution of the problem
end

"""
`BoundaryCondition` - abstract type to include boundary conditions
"""
abstract type BoundaryCondition end

"""
`NeumannBC` - structure to house a Neumann boundary condition, i.e., a boundary condition for which

uₓ(a) = vₗ  or  uₓ(b) = vᵣ
"""
struct NeumannBC <: BoundaryCondition
    boundary_condition::Real
end

"""
`PeriodicBC` - structure to house a Periodic boundary condition, i.e., a boundary condition for which

u(a) = u(b), uₓ(a) = uₓ(b)
"""
struct PeriodicBC <: BoundaryCondition
end

"""
`TimeSteppingMethod` - abstract type to include time stepping methods
"""
abstract type TimeSteppingMethod end

"""
`BackwardEulerMethod` - structure for Backward Euler method
"""
struct BackwardEulerMethod <: TimeSteppingMethod
end

"""
`CrankNicolsonMethod` - structure for Crank-Nicolson method
"""
struct CrankNicolsonMethod <: TimeSteppingMethod
end