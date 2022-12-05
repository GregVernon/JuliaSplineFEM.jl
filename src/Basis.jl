module Basis

function affineMapping( domain, target_domain, x )
    A = [ [ 1.0, 1.0 ] [ domain[1], domain[2] ] ]
    c = A \ target_domain
    fx =  c[1] + c[2] * x
    return fx
end

function evalBernstein( degree, basis_idx, x, domain )
    ξ = affineMapping( domain, [ 0.0, 1.0 ], x )
    term1 = binomial( degree, basis_idx )
    term2 = ξ ^ basis_idx
    term3 = ( 1.0 - ξ ) ^ ( degree - basis_idx )
    basis_val = term1 * term2 * term3
    return basis_val
end

end