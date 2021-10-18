module AllenCahn

using SparseArrays
using LinearAlgebra

include("structures.jl")
export AllenCahnProblem1D, NeumannBC, PeriodicBC, BackwardEulerMethod, CrankNicolsonMethod

include("constructors.jl")
export AllenCahnProblem1D

include("assembly.jl")
export assemble_system!

include("solve.jl")
export solve

end #module