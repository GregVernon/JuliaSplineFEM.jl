using BenchmarkTools
using Test

import Statistics 

include( "../src/Basis.jl" )

@testset "Affine Mapping Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    # Unit Domain --> Biunit Domain
    @test Basis.affineMapping( unit_domain, biunit_domain, unit_domain[1] ) ≈ biunit_domain[1]
    @test Basis.affineMapping( unit_domain, biunit_domain, unit_domain[2] ) ≈ biunit_domain[2]
    @test Basis.affineMapping( unit_domain, biunit_domain, 0.5 ) ≈ 0.0
    # Biunit Domain --> Unit Domain
    @test Basis.affineMapping( biunit_domain, unit_domain, biunit_domain[1] ) ≈ unit_domain[1]
    @test Basis.affineMapping( biunit_domain, unit_domain, biunit_domain[2] ) ≈ unit_domain[2]
    @test Basis.affineMapping( biunit_domain, unit_domain, 0.0 ) ≈ 0.5
    # Mesh Element Domain --> Biunit Domain
    @test Basis.affineMapping( elem_domain, biunit_domain, elem_domain[1] ) ≈ biunit_domain[1]
    @test Basis.affineMapping( elem_domain, biunit_domain, elem_domain[2] ) ≈ biunit_domain[2]
    @test Basis.affineMapping( elem_domain, biunit_domain, 17.0 ) ≈ 0.0
    # Mesh Element Domain --> Unit Domain
    @test Basis.affineMapping( elem_domain, unit_domain, elem_domain[1] ) ≈ unit_domain[1]
    @test Basis.affineMapping( elem_domain, unit_domain, elem_domain[2] ) ≈ unit_domain[2]
    @test Basis.affineMapping( elem_domain, unit_domain, 17.0 ) ≈ 0.5
    # Biunit Domain --> Mesh Element Domain
    @test Basis.affineMapping( biunit_domain, elem_domain, biunit_domain[1] ) ≈ elem_domain[1]
    @test Basis.affineMapping( biunit_domain, elem_domain, biunit_domain[2] ) ≈ elem_domain[2]
    @test Basis.affineMapping( biunit_domain, elem_domain, 0.0 ) ≈ 17.0
    # Unit Domain --> Mesh Element Domain
    @test Basis.affineMapping( unit_domain, elem_domain, unit_domain[1] ) ≈ elem_domain[1]
    @test Basis.affineMapping( unit_domain, elem_domain, unit_domain[2] ) ≈ elem_domain[2]
    @test Basis.affineMapping( unit_domain, elem_domain, 0.5 ) ≈ 17.0
end

@testset "Affine Mapping Performance Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    # Unit Domain --> Biunit Domain
    @test Statistics.median( @benchmark Basis.affineMapping( $unit_domain, $biunit_domain, 0.5 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $biunit_domain, $unit_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $elem_domain, $biunit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $elem_domain, $unit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $biunit_domain, $elem_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $unit_domain, $elem_domain, 0.5 ) ).time <= 50.0
end

@testset "Bernstein Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test Basis.evalBernstein( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test Basis.evalBernstein( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test Basis.evalBernstein( 1, 1, unit_domain, 0.5 ) ≈ 0.5
        @test Basis.evalBernstein( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test Basis.evalBernstein( 1, 1, unit_domain, 1.0 ) ≈ 0.0
        @test Basis.evalBernstein( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test Basis.evalBernstein( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test Basis.evalBernstein( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test Basis.evalBernstein( 1, 1, biunit_domain, +0.0 ) ≈ 0.5
        @test Basis.evalBernstein( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test Basis.evalBernstein( 1, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test Basis.evalBernstein( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test Basis.evalBernstein( 2, 1, unit_domain , 0.0 ) ≈ 1.00
        @test Basis.evalBernstein( 2, 2, unit_domain , 0.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 3, unit_domain , 0.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 1, unit_domain , 0.5 ) ≈ 0.25
        @test Basis.evalBernstein( 2, 2, unit_domain , 0.5 ) ≈ 0.50
        @test Basis.evalBernstein( 2, 3, unit_domain , 0.5 ) ≈ 0.25
        @test Basis.evalBernstein( 2, 1, unit_domain , 1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 2, unit_domain , 1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 3, unit_domain , 1.0 ) ≈ 1.00
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test Basis.evalBernstein( 2, 1, biunit_domain, -1.0 ) ≈ 1.00
        @test Basis.evalBernstein( 2, 2, biunit_domain, -1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 3, biunit_domain, -1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 1, biunit_domain, +0.0 ) ≈ 0.25
        @test Basis.evalBernstein( 2, 2, biunit_domain, +0.0 ) ≈ 0.50
        @test Basis.evalBernstein( 2, 3, biunit_domain, +0.0 ) ≈ 0.25
        @test Basis.evalBernstein( 2, 1, biunit_domain, +1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 2, biunit_domain, +1.0 ) ≈ 0.00
        @test Basis.evalBernstein( 2, 3, biunit_domain, +1.0 ) ≈ 1.00
    end
end

@testset "Lagrange Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test Basis.evalLagrange( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test Basis.evalLagrange( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test Basis.evalLagrange( 1, 1, unit_domain, 0.5 ) ≈ 0.5
        @test Basis.evalLagrange( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test Basis.evalLagrange( 1, 1, unit_domain, 1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test Basis.evalLagrange( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test Basis.evalLagrange( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 1, 1, biunit_domain, +0.0 ) ≈ 0.5
        @test Basis.evalLagrange( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test Basis.evalLagrange( 1, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test Basis.evalLagrange( 2, 1, unit_domain, 0.0 ) ≈ 1.0
        @test Basis.evalLagrange( 2, 2, unit_domain, 0.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 3, unit_domain, 0.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 1, unit_domain, 0.5 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 2, unit_domain, 0.5 ) ≈ 1.0
        @test Basis.evalLagrange( 2, 3, unit_domain, 0.5 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 1, unit_domain, 1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 2, unit_domain, 1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 3, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test Basis.evalLagrange( 2, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test Basis.evalLagrange( 2, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 3, biunit_domain, -1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 1, biunit_domain, +0.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 2, biunit_domain, +0.0 ) ≈ 1.0
        @test Basis.evalLagrange( 2, 3, biunit_domain, +0.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 1, biunit_domain, +1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 2, biunit_domain, +1.0 ) ≈ 0.0
        @test Basis.evalLagrange( 2, 3, biunit_domain, +1.0 ) ≈ 1.0
    end
end

@testset "Legendre Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test Basis.evalLegendre( 1, 1, unit_domain, 0.0 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 1, unit_domain, 0.5 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 1, unit_domain, 1.0 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 2, unit_domain, 0.0 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 0.0 )
        @test Basis.evalLegendre( 1, 2, unit_domain, 0.5 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 0.5 )
        @test Basis.evalLegendre( 1, 2, unit_domain, 1.0 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 1.0 )
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test Basis.evalLegendre( 1, 1, biunit_domain, -1.0 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 1, biunit_domain, +0.0 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 1, biunit_domain, +1.0 ) ≈ +1.0
        @test Basis.evalLegendre( 1, 2, biunit_domain, -1.0 ) ≈ -1.0
        @test Basis.evalLegendre( 1, 2, biunit_domain, +0.0 ) ≈ +0.0
        @test Basis.evalLegendre( 1, 2, biunit_domain, +1.0 ) ≈ +1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test Basis.evalLegendre( 2, 1, unit_domain, 0.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 1, unit_domain, 0.5 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 1, unit_domain, 1.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 2, unit_domain, 0.0 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 0.0 )
        @test Basis.evalLegendre( 2, 2, unit_domain, 0.5 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 0.5 )
        @test Basis.evalLegendre( 2, 2, unit_domain, 1.0 ) ≈ Basis.affineMapping( unit_domain, biunit_domain, 1.0 )
        @test Basis.evalLegendre( 2, 3, unit_domain, 0.0 ) ≈ 0.5 * ( 3.0 * Basis.affineMapping( unit_domain, biunit_domain, 0.0 ) ^ 2.0 - 1.0 )
        @test Basis.evalLegendre( 2, 3, unit_domain, 0.5 ) ≈ 0.5 * ( 3.0 * Basis.affineMapping( unit_domain, biunit_domain, 0.5 ) ^ 2.0 - 1.0 )
        @test Basis.evalLegendre( 2, 3, unit_domain, 1.0 ) ≈ 0.5 * ( 3.0 * Basis.affineMapping( unit_domain, biunit_domain, 1.0 ) ^ 2.0 - 1.0 )
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test Basis.evalLegendre( 2, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 1, biunit_domain, +0.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 1, biunit_domain, +1.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 2, biunit_domain, -1.0 ) ≈ -1.0
        @test Basis.evalLegendre( 2, 2, biunit_domain, +0.0 ) ≈ 0.0
        @test Basis.evalLegendre( 2, 2, biunit_domain, +1.0 ) ≈ 1.0
        @test Basis.evalLegendre( 2, 3, biunit_domain, -1.0 ) ≈ 0.5 * ( 3.0 * (-1.0)^2 - 1.0 )
        @test Basis.evalLegendre( 2, 3, biunit_domain, +0.0 ) ≈ 0.5 * ( 3.0 * (+0.0)^2 - 1.0 )
        @test Basis.evalLegendre( 2, 3, biunit_domain, +1.0 ) ≈ 0.5 * ( 3.0 * (+1.0)^2 - 1.0 )
    end
end

@testset "Monomial Basis Unit Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    @testset "Linear Basis, Unit Domain" begin
        @test Basis.evalMonomial( 1, 1, unit_domain, 0.0 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, unit_domain, 0.0 ) ≈ 0.0
        @test Basis.evalMonomial( 1, 1, unit_domain, 0.5 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, unit_domain, 0.5 ) ≈ 0.5
        @test Basis.evalMonomial( 1, 1, unit_domain, 1.0 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, unit_domain, 1.0 ) ≈ 1.0
    end
    @testset "Linear Basis, Biunit Domain" begin
        @test Basis.evalMonomial( 1, 1, biunit_domain, -1.0 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, biunit_domain, -1.0 ) ≈ 0.0
        @test Basis.evalMonomial( 1, 1, biunit_domain, +0.0 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, biunit_domain, +0.0 ) ≈ 0.5
        @test Basis.evalMonomial( 1, 1, biunit_domain, +1.0 ) ≈ 1.0
        @test Basis.evalMonomial( 1, 2, biunit_domain, +1.0 ) ≈ 1.0
    end
    @testset "Quadratic Basis, Unit Domain" begin
        @test Basis.evalMonomial( 2, 1, unit_domain, 0.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, unit_domain, 0.0 ) ≈ 0.00
        @test Basis.evalMonomial( 2, 3, unit_domain, 0.0 ) ≈ 0.00
        @test Basis.evalMonomial( 2, 1, unit_domain, 0.5 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, unit_domain, 0.5 ) ≈ 0.50
        @test Basis.evalMonomial( 2, 3, unit_domain, 0.5 ) ≈ 0.25
        @test Basis.evalMonomial( 2, 1, unit_domain, 1.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, unit_domain, 1.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 3, unit_domain, 1.0 ) ≈ 1.00
    end
    @testset "Quadratic Basis, Biunit Domain" begin
        @test Basis.evalMonomial( 2, 1, biunit_domain, -1.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, biunit_domain, -1.0 ) ≈ 0.00
        @test Basis.evalMonomial( 2, 3, biunit_domain, -1.0 ) ≈ 0.00
        @test Basis.evalMonomial( 2, 1, biunit_domain, +0.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, biunit_domain, +0.0 ) ≈ 0.50
        @test Basis.evalMonomial( 2, 3, biunit_domain, +0.0 ) ≈ 0.25
        @test Basis.evalMonomial( 2, 1, biunit_domain, +1.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 2, biunit_domain, +1.0 ) ≈ 1.00
        @test Basis.evalMonomial( 2, 3, biunit_domain, +1.0 ) ≈ 1.00
    end
end

@testset "Compute Legendre Roots" begin
    @test Basis.computeLegendreRoots( 0 ) == empty(Vector{Float64}([]))
    @test isapprox( Basis.computeLegendreRoots( 1 ), 0.0 )
    @test isapprox( Basis.computeLegendreRoots( 2 ), [ -1/sqrt(3), 1/sqrt(3) ] )
    @test isapprox( Basis.computeLegendreRoots( 3 ), [ -sqrt(3/5), 0, sqrt(3/5) ] )
end