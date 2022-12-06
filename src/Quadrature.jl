module Quadrature
import ..Basis
import LinearAlgebra 

function integrateOverElement( fun, domain, method )
end

function getGaussLegendreQuadrature( num_points )
    ξ = Basis.computeLegendreRoots( num_points )
    w = solveLinearMomentFit( ξ )
    return ξ, w
end

function solveLinearMomentFit( ξ )
    M, A = assembleLinearMomentFit( ξ )
    w = A \ M
    return w
end

function assembleLinearMomentFit( ξ )
    num_points = length( ξ )
    M = [ 2.0; zeros( 2 * num_points - 1 ) ]
    A = zeros( length( M ), num_points )
    for p = 1:length( M )
        for i = 1:num_points
            A[p, i] = Basis.evalLegendre( p, p, [ -1.0, 1.0 ], ξ[i] )
        end
    end
    return M, A
end

# struct quadratureParams()
#     method_name :: string
#     integration_domain :: Vector
#     integration_order :: Int
# end

end