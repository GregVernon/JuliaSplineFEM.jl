@testset "Affine Mapping Performance Test" begin
    unit_domain = [ 0.0, 1.0 ]
    biunit_domain = [ -1.0, 1.0 ]
    elem_domain = [ 15.5, 18.5 ]
    # Unit Domain --> Biunit Domain
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $unit_domain, $biunit_domain, 0.5 ) ).time <= 50.0
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $biunit_domain, $unit_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $elem_domain, $biunit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $elem_domain, $unit_domain, 17.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $biunit_domain, $elem_domain, 0.0 ) ).time <= 50.0
    @test Statistics.median( @benchmark JuliaSplineFEM.affineMapping( $unit_domain, $elem_domain, 0.5 ) ).time <= 50.0
end