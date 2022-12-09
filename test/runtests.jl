using BenchmarkTools
using Test
import Statistics 

include( "../src/JuliaSplineFEM.jl" )

@testset "Basis" begin
    include( "test_Basis.jl" )
end

@testset "Quadrature" begin
    include( "test_Quadrature.jl" )
end