struct AllenCahnProblem1D{TLBC, TRBC, TA}
    a::Real # Left endpoint
    b::Real # Right endpoint
    ϵ::Real # Model parameter
    nₓ::Int # Number of spatial mesh points
    nₜ::Int # Number of timesteps
    Δt::Real # Length of timestep
    left_bc::TLBC # Boundary condition on the left endpoint
    right_bc::TRBC # Boundary condition on the right endpoint
    A::TA # Matrix approximator of differential operator for the problem
    rhs::Vector # Right hand side of the end system
end


struct NeumannBC
    boundary_condition::Real
end