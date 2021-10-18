using AllenCahn
using Test
using LinearAlgebra

@testset "Neumann" begin
    @test include("neumanntest1.jl")
    @test include("neumanntest2.jl")
end

@testset "Periodic" begin
    @test include("periodictest1.jl")
    @test include("periodictest2.jl")
end