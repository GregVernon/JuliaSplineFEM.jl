@testset "Affine Mapping Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    # Unit Domain --> Biunit Domain
    @test JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, unit_domain[1] ) ≈ biunit_domain[1]
    @test JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, unit_domain[2] ) ≈ biunit_domain[2]
    @test JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.5 ) ≈ 0.0
    # Biunit Domain --> Unit Domain
    @test JuliaSplineFEM.affineMapping( biunit_domain, unit_domain, biunit_domain[1] ) ≈ unit_domain[1]
    @test JuliaSplineFEM.affineMapping( biunit_domain, unit_domain, biunit_domain[2] ) ≈ unit_domain[2]
    @test JuliaSplineFEM.affineMapping( biunit_domain, unit_domain, 0.0 ) ≈ 0.5
    # Mesh Element Domain --> Biunit Domain
    @test JuliaSplineFEM.affineMapping( elem_domain, biunit_domain, elem_domain[1] ) ≈ biunit_domain[1]
    @test JuliaSplineFEM.affineMapping( elem_domain, biunit_domain, elem_domain[2] ) ≈ biunit_domain[2]
    @test JuliaSplineFEM.affineMapping( elem_domain, biunit_domain, 17.0 ) ≈ 0.0
    # Mesh Element Domain --> Unit Domain
    @test JuliaSplineFEM.affineMapping( elem_domain, unit_domain, elem_domain[1] ) ≈ unit_domain[1]
    @test JuliaSplineFEM.affineMapping( elem_domain, unit_domain, elem_domain[2] ) ≈ unit_domain[2]
    @test JuliaSplineFEM.affineMapping( elem_domain, unit_domain, 17.0 ) ≈ 0.5
    # Biunit Domain --> Mesh Element Domain
    @test JuliaSplineFEM.affineMapping( biunit_domain, elem_domain, biunit_domain[1] ) ≈ elem_domain[1]
    @test JuliaSplineFEM.affineMapping( biunit_domain, elem_domain, biunit_domain[2] ) ≈ elem_domain[2]
    @test JuliaSplineFEM.affineMapping( biunit_domain, elem_domain, 0.0 ) ≈ 17.0
    # Unit Domain --> Mesh Element Domain
    @test JuliaSplineFEM.affineMapping( unit_domain, elem_domain, unit_domain[1] ) ≈ elem_domain[1]
    @test JuliaSplineFEM.affineMapping( unit_domain, elem_domain, unit_domain[2] ) ≈ elem_domain[2]
    @test JuliaSplineFEM.affineMapping( unit_domain, elem_domain, 0.5 ) ≈ 17.0
end

@testset "Bernstein Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalBernstein( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalBernstein( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalBernstein( 1, 1, unit_domain, 0.5 ) ≈ 0.5
        @test JuliaSplineFEM.evalBernstein( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test JuliaSplineFEM.evalBernstein( 1, 1, unit_domain, 1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalBernstein( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalBernstein( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalBernstein( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalBernstein( 1, 1, biunit_domain, +0.0 ) ≈ 0.5
        @test JuliaSplineFEM.evalBernstein( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test JuliaSplineFEM.evalBernstein( 1, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalBernstein( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalBernstein( 2, 1, unit_domain , 0.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalBernstein( 2, 2, unit_domain , 0.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 3, unit_domain , 0.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 1, unit_domain , 0.5 ) ≈ 0.25
        @test JuliaSplineFEM.evalBernstein( 2, 2, unit_domain , 0.5 ) ≈ 0.50
        @test JuliaSplineFEM.evalBernstein( 2, 3, unit_domain , 0.5 ) ≈ 0.25
        @test JuliaSplineFEM.evalBernstein( 2, 1, unit_domain , 1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 2, unit_domain , 1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 3, unit_domain , 1.0 ) ≈ 1.00
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalBernstein( 2, 1, biunit_domain, -1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalBernstein( 2, 2, biunit_domain, -1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 3, biunit_domain, -1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 1, biunit_domain, +0.0 ) ≈ 0.25
        @test JuliaSplineFEM.evalBernstein( 2, 2, biunit_domain, +0.0 ) ≈ 0.50
        @test JuliaSplineFEM.evalBernstein( 2, 3, biunit_domain, +0.0 ) ≈ 0.25
        @test JuliaSplineFEM.evalBernstein( 2, 1, biunit_domain, +1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 2, biunit_domain, +1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalBernstein( 2, 3, biunit_domain, +1.0 ) ≈ 1.00
    end
end

@testset "Lagrange Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalLagrange( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 1, 1, unit_domain, 0.5 ) ≈ 0.5
        @test JuliaSplineFEM.evalLagrange( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test JuliaSplineFEM.evalLagrange( 1, 1, unit_domain, 1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalLagrange( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 1, 1, biunit_domain, +0.0 ) ≈ 0.5
        @test JuliaSplineFEM.evalLagrange( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test JuliaSplineFEM.evalLagrange( 1, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalLagrange( 2, 1, unit_domain, 0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, unit_domain, 0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, unit_domain, 0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 1, unit_domain, 0.5 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, unit_domain, 0.5 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, unit_domain, 0.5 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 1, unit_domain, 1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, unit_domain, 1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalLagrange( 2, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, biunit_domain, -1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 1, biunit_domain, +0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, biunit_domain, +0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, biunit_domain, +0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 2, biunit_domain, +1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLagrange( 2, 3, biunit_domain, +1.0 ) ≈ 1.0
    end
end

@testset "Legendre Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalLegendre( 1, 1, unit_domain, 0.0 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 1, unit_domain, 0.5 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 1, unit_domain, 1.0 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 2, unit_domain, 0.0 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.0 )
        @test JuliaSplineFEM.evalLegendre( 1, 2, unit_domain, 0.5 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.5 )
        @test JuliaSplineFEM.evalLegendre( 1, 2, unit_domain, 1.0 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 1.0 )
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalLegendre( 1, 1, biunit_domain, -1.0 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 1, biunit_domain, +0.0 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 1, biunit_domain, +1.0 ) ≈ +1.0
        @test JuliaSplineFEM.evalLegendre( 1, 2, biunit_domain, -1.0 ) ≈ -1.0
        @test JuliaSplineFEM.evalLegendre( 1, 2, biunit_domain, +0.0 ) ≈ +0.0
        @test JuliaSplineFEM.evalLegendre( 1, 2, biunit_domain, +1.0 ) ≈ +1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalLegendre( 2, 1, unit_domain, 0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 1, unit_domain, 0.5 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 1, unit_domain, 1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 2, unit_domain, 0.0 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 2, unit_domain, 0.5 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.5 )
        @test JuliaSplineFEM.evalLegendre( 2, 2, unit_domain, 1.0 ) ≈ JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 1.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 3, unit_domain, 0.0 ) ≈ 0.5 * ( 3.0 * JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.0 ) ^ 2.0 - 1.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 3, unit_domain, 0.5 ) ≈ 0.5 * ( 3.0 * JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 0.5 ) ^ 2.0 - 1.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 3, unit_domain, 1.0 ) ≈ 0.5 * ( 3.0 * JuliaSplineFEM.affineMapping( unit_domain, biunit_domain, 1.0 ) ^ 2.0 - 1.0 )
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalLegendre( 2, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 1, biunit_domain, +0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 1, biunit_domain, +1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 2, biunit_domain, -1.0 ) ≈ -1.0
        @test JuliaSplineFEM.evalLegendre( 2, 2, biunit_domain, +0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalLegendre( 2, 2, biunit_domain, +1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalLegendre( 2, 3, biunit_domain, -1.0 ) ≈ 0.5 * ( 3.0 * (-1.0)^2 - 1.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 3, biunit_domain, +0.0 ) ≈ 0.5 * ( 3.0 * (+0.0)^2 - 1.0 )
        @test JuliaSplineFEM.evalLegendre( 2, 3, biunit_domain, +1.0 ) ≈ 0.5 * ( 3.0 * (+1.0)^2 - 1.0 )
    end
end

@testset "Monomial Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalMonomial( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalMonomial( 1, 1, unit_domain, 0.5 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test JuliaSplineFEM.evalMonomial( 1, 1, unit_domain, 1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalMonomial( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test JuliaSplineFEM.evalMonomial( 1, 1, biunit_domain, +0.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test JuliaSplineFEM.evalMonomial( 1, 1, biunit_domain, +1.0 ) ≈ 1.0
        @test JuliaSplineFEM.evalMonomial( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test JuliaSplineFEM.evalMonomial( 2, 1, unit_domain, 0.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, unit_domain, 0.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalMonomial( 2, 3, unit_domain, 0.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalMonomial( 2, 1, unit_domain, 0.5 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, unit_domain, 0.5 ) ≈ 0.50
        @test JuliaSplineFEM.evalMonomial( 2, 3, unit_domain, 0.5 ) ≈ 0.25
        @test JuliaSplineFEM.evalMonomial( 2, 1, unit_domain, 1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, unit_domain, 1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 3, unit_domain, 1.0 ) ≈ 1.00
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test JuliaSplineFEM.evalMonomial( 2, 1, biunit_domain, -1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, biunit_domain, -1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalMonomial( 2, 3, biunit_domain, -1.0 ) ≈ 0.00
        @test JuliaSplineFEM.evalMonomial( 2, 1, biunit_domain, +0.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, biunit_domain, +0.0 ) ≈ 0.50
        @test JuliaSplineFEM.evalMonomial( 2, 3, biunit_domain, +0.0 ) ≈ 0.25
        @test JuliaSplineFEM.evalMonomial( 2, 1, biunit_domain, +1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 2, biunit_domain, +1.0 ) ≈ 1.00
        @test JuliaSplineFEM.evalMonomial( 2, 3, biunit_domain, +1.0 ) ≈ 1.00
    end
end

@testset "Compute Legendre Roots" begin
    @test JuliaSplineFEM.computeLegendreRoots( 0 ) == empty(Vector{Float64}([]))
    @test isapprox( JuliaSplineFEM.computeLegendreRoots( 1 ), [ 0.0 ] )
    @test isapprox( JuliaSplineFEM.computeLegendreRoots( 2 ), [ -1/sqrt(3), 1/sqrt(3) ] )
    @test isapprox( JuliaSplineFEM.computeLegendreRoots( 3 ), [ -sqrt(3/5), 0, sqrt(3/5) ] )
end

@testset "polynomialChangeOfBasis" begin
    for n = 0:10
        domain = sort( rand(2)*2 .- 1 )
        c = rand( n+1 )*2 .- 1
        d, _ = JuliaSplineFEM.polynomialChangeOfBasis( JuliaSplineFEM.evalLegendre, JuliaSplineFEM.evalBernstein, c, domain )
        pc = (x) -> transpose( c ) * JuliaSplineFEM.evalElementBasis( JuliaSplineFEM.evalLegendre, n, domain, x )
        pd = (x) -> transpose( d ) * JuliaSplineFEM.evalElementBasis( JuliaSplineFEM.evalBernstein, n, domain, x )
        num_qp = JuliaSplineFEM.computeNumGaussPointsFromPolynomialDegree( n ) + 1
        @test isapprox( JuliaSplineFEM.integrateOverElement( (x)->abs( pc(x) - pd(x) ), domain, num_qp ), 0.0, atol = 1e-12 )
    end
end