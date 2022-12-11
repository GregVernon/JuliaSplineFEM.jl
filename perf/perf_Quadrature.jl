@testset "Perf: integrateOverElement" begin
    biunit_domain = [ -1.0, 1.0 ]
    fun = (x) -> cos(x)
    println( Statistics.median( @benchmark JuliaSplineFEM.integrateOverElement( $fun, $biunit_domain, 4 ) ).time )
end