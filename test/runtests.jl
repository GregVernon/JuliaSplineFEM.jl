using JuliaSplineFEM
using Test

@testset "JuliaSplineFEM.jl" begin
    include( "test_Basis.jl" )
    include( "test_Quadrature.jl" )
end