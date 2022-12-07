using BenchmarkTools
using Test

import Statistics 

include( "../src/Quadrature.jl" )
# import ..JuliaSplineFEM
import .Quadrature

@testset "Integrate Polynomials Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    
    for degree = 0:10
        num_qp = Quadrature.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = 1.0 / ( degree + 1 )
        @test isapprox( Quadrature.integrateOverElement( (x) -> x^degree, unit_domain, num_qp ), exact_int, atol = 1e-12 )
    end
    
    for degree = 0:10
        num_qp = Quadrature.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = ( ( 1.0 ^ ( degree + 1 ) ) - ( (-1.0) ^ ( degree + 1 ) ) ) / ( degree + 1 )
        @test isapprox( Quadrature.integrateOverElement( (x) -> x^degree, biunit_domain, num_qp ), exact_int, atol = 1e-12 )
    end

    for degree = 0:10
        num_qp = Quadrature.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = ( ( elem_domain[2] ^ ( degree + 1 ) ) - ( elem_domain[1] ^ ( degree + 1 ) ) ) / ( degree + 1 )
        @test isapprox( Quadrature.integrateOverElement( (x) -> x^degree, elem_domain, num_qp ), exact_int, atol = 1e-12, rtol = 1e-12 )
    end
    
end