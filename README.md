# AllenCahn.jl
Numerical Solver Package for the one-dimensional deterministic Allen-Cahn equation

This package solves the one-dimensional deterministic Allen-Cahn equation $\newline \newline$
$$
u_{t} = \epsilon u_{xx} + u - u^{3} \;\; \text{on } (a, b)
\newline
$$
where $\epsilon > 0$ is user-specified. The package supports use of (constant) Neumann and Periodic boundary conditions, and implements the semi-implicit Euler and Crank-Nicolson time stepping schemes. The spatial discretization scheme is simple finite differences using sparse matrices.
