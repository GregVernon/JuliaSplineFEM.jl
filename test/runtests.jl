using BenchmarkTools
using Test
import Statistics 

@testset "JuliaSplineFEM.jl" begin
    include( "../src/JuliaSplineFEM.jl" )
    include( "test_Basis.jl" )
    include( "test_Quadrature.jl" )
end