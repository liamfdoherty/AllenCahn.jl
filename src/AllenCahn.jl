module AllenCahn

using SparseArrays
using LinearAlgebra

include("structures.jl")
export AllenCahnProblem1D, NeumannBC

include("constructors.jl")
export AllenCahnProblem1D

end #module