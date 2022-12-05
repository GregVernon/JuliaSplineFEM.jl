module Basis

function affineMapping( domain, target_domain, x )
    x -= domain[1]
    x *= ( target_domain[2] - target_domain[1] ) / ( domain[2] - domain[1] )
    x += target_domain[1]
    return x
end

function evalBernstein( degree, basis_idx, domain, x )
    ξ = affineMapping( domain, [ 0.0, 1.0 ], x )
    term1 = binomial( degree, basis_idx - 1 )
    term2 = ξ ^ ( basis_idx - 1 )
    term3 = ( 1.0 - ξ ) ^ ( degree - ( basis_idx - 1 ) )
    basis_val = term1 * term2 * term3
    return basis_val
end

end