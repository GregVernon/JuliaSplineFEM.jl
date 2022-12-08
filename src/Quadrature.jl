module Quadrature

import LinearAlgebra 
import ForwardDiff

push!(LOAD_PATH, "./src" )
import Basis

function integrateOverElement( fun, domain, num_points )
    ξ, w = getGaussLegendreQuadrature( num_points )
    I = 0.0
    for i = axes( ξ, 1 )
        x = Basis.affineMapping( [-1, 1], domain, ξ[i] )
        J = ForwardDiff.derivative( (ξ) -> Basis.affineMapping( [-1, 1], domain, ξ ) , ξ[i] )
        I += fun( x ) * w[i] * J
    end
    return I
end

function getGaussLegendreQuadrature( num_points )
    ξ = Basis.computeLegendreRoots( num_points )
    ξ = isempty( ξ ) ? [0.0] : ξ
    w = solveLinearMomentFit( ξ )
    return ξ, w
end

function computeNumGaussPointsFromPolynomialDegree( degree )
    return Int( ceil( ( degree + 1.0 ) / 2.0 ) )
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