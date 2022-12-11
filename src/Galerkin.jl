function computeGalerkinApproximation( uspline::USpline, f )
    M = assembleGramMatrix( uspline )
    F = assembleForceVector( uspline, f )
    d = M \ F
    return d
end

function assembleGramMatrix( uspline::USpline )
    num_nodes = getNumNodes( uspline )
    num_elems = getNumElems( uspline )
    gram_matrix = zeros( num_nodes, num_nodes )
    for elem_idx = 1 : num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_domain = getElementDomain( uspline, elem_id )
        elem_degree = getElementDegree( uspline, elem_id )
        elem_nodes = getElementNodeIds( uspline, elem_id )
        elem_extraction_operator = getElementExtractionOperator( uspline, elem_id )
        num_qp = computeNumGaussPointsFromPolynomialDegree( elem_degree * 2 )
        for local_node_id_1 = 1 : length( elem_nodes )
            global_node_id_1 = elem_nodes[ local_node_id_1 ]
            N1 = (x) -> evalSplineBasis( elem_extraction_operator, local_node_id_1, elem_domain, x )
            for local_node_id_2 = 1 : length( elem_nodes )
                global_node_id_2 = elem_nodes[ local_node_id_2 ]
                N2 = (x) -> evalSplineBasis( elem_extraction_operator, local_node_id_2, elem_domain, x )
                integrand = (x) -> N1(x) * N2(x)
                gram_matrix[ global_node_id_1, global_node_id_2 ] += integrateOverElement( integrand, elem_domain, num_qp )
            end
        end
    end
    return gram_matrix
end

function assembleForceVector( uspline::USpline, f )
    num_nodes = getNumNodes( uspline )
    num_elems = getNumElems( uspline )
    force_vector = zeros( num_nodes )
    for elem_idx = 1 : num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_domain = getElementDomain( uspline, elem_id )
        elem_degree = getElementDegree( uspline, elem_id )
        elem_nodes = getElementNodeIds( uspline, elem_id )
        elem_extraction_operator = getElementExtractionOperator( uspline, elem_id )
        num_qp = computeNumGaussPointsFromPolynomialDegree( elem_degree * 2 )
        for local_node_id_1 = 1 : length( elem_nodes )
            global_node_id_1 = elem_nodes[ local_node_id_1 ]
            N1 = (x) -> evalSplineBasis( elem_extraction_operator, local_node_id_1, elem_domain, x )
            integrand = (x) -> N1(x) * f( x )
            force_vector[ global_node_id_1 ] += integrateOverElement( integrand, elem_domain, num_qp )
        end
    end
    return force_vector
end

function evaluateSolutionAt( uspline::USpline, sol_coeff, x )
    elem_id = getElementIdContainingPoint( uspline, x )
    elem_nodes = getElementNodeIds( uspline, elem_id )
    elem_domain = getElementDomain( uspline, elem_id )
    elem_degree = getElementDegree( uspline, elem_id )
    elem_extraction_operator = getElementExtractionOperator( uspline, elem_id )
    sol = 0.0
    for local_node_id_1 = 1 : length( elem_nodes )
        global_node_id_1 = elem_nodes[ local_node_id_1 ]
        sol += sol_coeff[ global_node_id_1 ] * evalSplineBasis( elem_extraction_operator, local_node_id_1, elem_domain, x )
    end
    return sol
end