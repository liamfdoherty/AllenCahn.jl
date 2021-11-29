# AllenCahn.jl
Numerical Solver Package for the one-dimensional deterministic Allen-Cahn equation

This package solves the one-dimensional deterministic Allen-Cahn equation

![equation](https://latex.codecogs.com/png.latex?u_%7Bt%7D%20%3D%20%5Cepsilon%20u_%7Bxx%7D%20&plus;%20u%20-%20u%5E%7B3%7D%20%5C%3B%5C%3B%20%5Ctext%7Bon%20%7D%20%28a%2C%20b%29)

where ![equation](https://latex.codecogs.com/png.latex?%5Cbg_white%20%5Cepsilon%20%3E%200) is user-specified. The package supports use of (constant) Neumann and Periodic boundary conditions, and implements the semi-implicit Euler and Crank-Nicolson time stepping schemes. The spatial discretization scheme is simple finite differences using sparse matrices.
