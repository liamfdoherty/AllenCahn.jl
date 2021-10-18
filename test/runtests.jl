using AllenCahn
using Test
using LinearAlgebra

@testset "Neumann" begin
    @test include("neumanntest1.jl")
end

@testset "Periodic" begin
    @test include("periodictest1.jl")
end