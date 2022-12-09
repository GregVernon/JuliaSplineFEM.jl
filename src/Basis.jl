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

function evalLagrange( degree, basis_idx, domain, x )
    ξ = affineMapping( domain, [ -1.0, 1.0 ], x )
    nodes = LinRange( -1.0, 1.0, degree + 1 )
    basis_val = 1.0
    for i = 1 : degree + 1
        if i != basis_idx
            basis_val *= ( ξ - nodes[i] ) / ( nodes[basis_idx] - nodes[i] )
        end
    end
    return basis_val
end

function evalLegendre( degree, basis_idx, domain, x )
    ξ = affineMapping( domain, [ -1.0, 1.0 ], x )
    if basis_idx == 1
        basis_val = 1.0
    elseif basis_idx == 2
        basis_val = ξ
    else
        n = basis_idx - 2
        term1 = ( ( 2 * n ) + 1 ) * ξ * evalLegendre( n, n + 1, domain, x )
        term2 = n * evalLegendre( n-1, n, domain, x )
        basis_val = ( term1 - term2 ) / ( n + 1 )
    end
    return basis_val
end

function evalMonomial( degree, basis_idx, domain, x )
    ξ = affineMapping( domain, [ 0.0, 1.0 ], x )
    return ξ ^ ( basis_idx - 1 )
end

function evalElementBasis( basis::Function, degree, domain, x )
    f = zeros( degree + 1 )
    for i = 1:degree+1
        f[i] = basis( degree, i, domain, x )
    end
    return f
end

function polynomialChangeOfBasis( poly_basis_1, poly_basis_2, coeff_1, domain )
    degree = length( coeff_1 ) - 1
    num_qp = computeNumGaussPointsFromPolynomialDegree( 2 * degree )
    C = zeros( degree + 1, degree + 1 )
    D = zeros( degree + 1, degree + 1 )
    for i = 1 : degree + 1
        N2i = (x) -> poly_basis_2( degree, i, domain, x )
        for j = 1 : degree + 1
            N2j = (x) -> poly_basis_2( degree, j, domain, x )
            N1j = (x) -> poly_basis_1( degree, j, domain, x )
            C[i,j] += integrateOverElement( (x)-> N2i(x) * N1j(x), domain, num_qp )
            D[i,j] += integrateOverElement( (x)-> N2i(x) * N2j(x), domain, num_qp )
        end
    end

    R = D \ C
    d = R * coeff_1
    return d, R
end

function computeLegendreRoots( degree )
    if degree == 0
        roots = empty(Vector{Float64}([]))
    else
        nodes = LinRange( -1.0, 1.0, degree + 1 )
        V = SpecialMatrices.Vandermonde( nodes )
        b = zeros( Float64, degree + 1 )
        for i = axes( nodes, 1 )
            b[i] = evalLegendre( degree, degree + 1, [ -1.0, 1.0 ], nodes[i] )
        end
        coeff = V \ b
        P = Polynomials.Polynomial( coeff )
        C = SpecialMatrices.Companion( P )
        roots = LinearAlgebra.eigvals( C )
    end
    return roots
end