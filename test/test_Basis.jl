using BenchmarkTools
using Test

import Statistics 

include( "../src/Basis.jl" )

@testset "Affine Mapping Unit Test" begin
    unit_domain = [ 0, 1 ]
    biunit_domain = [ -1, 1 ]
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
    unit_domain = [ 0, 1 ]
    biunit_domain = [ -1, 1 ]
    elem_domain = [ 15.5, 18.5 ]
    # Unit Domain --> Biunit Domain
    @test Statistics.median( @benchmark Basis.affineMapping( $unit_domain, $biunit_domain, 0.5 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $biunit_domain, $unit_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $elem_domain, $biunit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $elem_domain, $unit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $biunit_domain, $elem_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark Basis.affineMapping( $unit_domain, $elem_domain, 0.5 ) ).time <= 50.0
end

