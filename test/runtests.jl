using AllenCahn
using Test
using LinearAlgebra

@testset "Neumann" begin
    @test include("neumanntest1.jl")
end