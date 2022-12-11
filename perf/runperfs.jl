using BenchmarkTools
using Test
import Statistics 

include( "../src/JuliaSplineFEM.jl" )

@testset "Basis" begin
    include( "perf_Basis.jl" )
end

@testset "Quadrature" begin
    include( "perf_Quadrature.jl" )
end

