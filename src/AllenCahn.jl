module AllenCahn

using SparseArrays
using LinearAlgebra

include("structures.jl")
export AllenCahnProblem1D, NeumannBC

include("constructors.jl")
export AllenCahnProblem1D

include("assembly.jl")
export assemble_system!

end #module