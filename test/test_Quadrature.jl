@testset "Integrate Polynomials Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    
    for degree = 0:10
        num_qp = JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = 1.0 / ( degree + 1 )
        @test isapprox( JuliaSplineFEM.integrateOverElement( (x) -> x^degree, unit_domain, num_qp ), exact_int, atol = 1e-12 )
    end
    
    for degree = 0:10
        num_qp = JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = ( ( 1.0 ^ ( degree + 1 ) ) - ( (-1.0) ^ ( degree + 1 ) ) ) / ( degree + 1 )
        @test isapprox( JuliaSplineFEM.integrateOverElement( (x) -> x^degree, biunit_domain, num_qp ), exact_int, atol = 1e-12 )
    end

    for degree = 0:10
        num_qp = JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( degree )
        exact_int = ( ( elem_domain[2] ^ ( degree + 1 ) ) - ( elem_domain[1] ^ ( degree + 1 ) ) ) / ( degree + 1 )
        @test isapprox( JuliaSplineFEM.integrateOverElement( (x) -> x^degree, elem_domain, num_qp ), exact_int, atol = 1e-12, rtol = 1e-12 )
    end
end

@testset "computeNumGaussPointsFromPolynomialDegree" begin
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 0 ) == 1
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 1 ) == 1
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 2 ) == 2
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 3 ) == 2
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 4 ) == 3
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 5 ) == 3
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 6 ) == 4
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 7 ) == 4
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 8 ) == 5
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 9 ) == 5
    @test JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( 10 ) == 6
end