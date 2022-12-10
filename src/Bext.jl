function readBEXT( filename )
    b = JSON.parsefile( filename )
    is_rational = b["is_rational"]
    num_nodes = b["num_nodes"]
    num_vertices = b["num_vertices"]
    num_elems = b["num_elems"]
    dim = b["dim"]
    num_coeff_vecs = b["num_coeff_vecs"]
    patch_id = b["patch_id"]
    elements = b["elements"]
    block_elem_map = []
    for i = axes( b["block_elem_map"], 1 )
        block_name = b["block_elem_map"][i]["block_name"]
        block_id = b["block_elem_map"][i]["block_id"]
        elem_ids = b["block_elem_map"][i]["elem_ids"]
        push!( block_elem_map, BlockElemMap( block_name, block_id, elem_ids ) )
    end
    vertex_connectivity = b["vertex_connectivity"]
    nodes = mapreduce(permutedims, vcat, b["nodes"])
    # Process Coefficients
    num_dense_coeff_blocks = b["coefficients"]["num_dense_coeff_blocks"]
    dense_coefficient_vectors = []
    for i = axes( b["coefficients"]["dense_coefficient_vectors"], 1 )
        push!( dense_coefficient_vectors, b["coefficients"]["dense_coefficient_vectors"][i]["components"] )
    end
    num_sparse_coeff_blocks = b["coefficients"]["num_sparse_coeff_blocks"]
    sparse_coefficient_vectors = []
    for i = axes( b["coefficients"]["sparse_coefficient_vectors"], 1 )
        push!( sparse_coefficient_vectors, b["coefficients"]["sparse_coefficient_vectors"][i]["components"] )
    end
    sparse_vector_info = b["coefficients"]["sparse_vector_info"]
    dense_vector_info = b["coefficients"]["dense_vector_info"]
    coefficients = Coefficients( num_dense_coeff_blocks, dense_coefficient_vectors, sparse_coefficient_vectors, num_sparse_coeff_blocks, sparse_vector_info, dense_vector_info )
    B = USpline( is_rational, num_nodes, num_vertices, num_elems, dim, num_coeff_vecs, patch_id, elements, block_elem_map, vertex_connectivity, nodes, coefficients)
    return B
end

struct BlockElemMap
    block_name::String
    block_id::Int
    elem_ids::Array{Int}
end

struct Coefficients
    num_dense_coeff_blocks::Int
    dense_coefficient_vectors::Array
    sparse_coefficient_vectors::Array
    num_sparse_coeff_blocks::Int
    sparse_vector_info::Array
    dense_vector_info::Array
end

struct USpline
    is_rational::Bool
    num_nodes::Int
    num_vertices::Int
    num_elems::Int
    dim::Int
    num_coeff_vecs::Int
    patch_id::Int
    elements::Dict{String, Any}
    block_elem_map::Vector{BlockElemMap}
    vertex_connectivity::Vector{Vector{Int}}
    nodes::Matrix
    coefficients::Coefficients
end

function getNumElems( uspline::USpline )
    return uspline.num_elems
end

function getNumVertices( uspline::USpline )
    return uspline.num_vertices
end

function getNumNodes( uspline::USpline )
    return size( getSplineNodes( uspline ), 1 )
end

function getNumBezierNodes( uspline::USpline )
    num_elems = getNumElems( uspline )
    num_bez_nodes = 0
    for elem_idx = 1:num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_degree = getElementDegree( uspline, elem_id )
        num_bez_nodes += elem_degree + 1
    end
    return num_bez_nodes
end

function getSplineNodes( uspline::USpline )
    return uspline.nodes
end

function getDomain( uspline::USpline )
    nodes = getSplineNodes( uspline )
    return [ minimum( nodes[:,1] ), maximum( nodes[:,1] ) ]
end

function elemIdFromElemIdx( uspline::USpline, elem_idx )
    element_blocks = uspline.elements["element_blocks"]
    elem_id = element_blocks[ elem_idx ]["us_cid"]
    return elem_id
end

function elemIdxFromElemId( uspline::USpline, elem_id )
    element_blocks = uspline.elements["element_blocks"]
    for elem_idx = axes( element_blocks, 1 )
        if element_blocks[ elem_idx ]["us_cid"] == elem_id
            return elem_idx
        end
    end
end

function getElementDegree( uspline::USpline, elem_id )
    elem_idx = elemIdxFromElemId( uspline, elem_id )
    return length( uspline.elements["element_blocks"][elem_idx]["node_ids"] ) - 1
end

function getElementDomain( uspline::USpline, elem_id )
    elem_bezier_nodes = getElementBezierNodes( uspline, elem_id )
    elem_domain = [ minimum( elem_bezier_nodes[:,1] ), maximum( elem_bezier_nodes[:,1] ) ]
    return elem_domain
end

function getElementNodeIds( uspline::USpline, elem_id )
    elem_idx = elemIdxFromElemId( uspline, elem_id )
    elem_node_ids = uspline.elements["element_blocks"][elem_idx]["node_ids"] .+ 1
    return elem_node_ids
end

function getElementBezierNodeIds( uspline::USpline, elem_id )
    num_elems = getNumElems( uspline )
    num_bez_nodes = 1
    for elem_idx = 1 : num_elems
        curr_elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_degree = getElementDegree( uspline, curr_elem_id )
        if elem_id == curr_elem_id
            elem_bez_node_ids = collect( num_bez_nodes : num_bez_nodes + elem_degree )
            return elem_bez_node_ids
        end
        num_bez_nodes += elem_degree + 1
    end
end

function getElementNodes( uspline::USpline, elem_id )
    elem_node_ids = getElementNodeIds( uspline, elem_id )
    spline_nodes = getSplineNodes( uspline )
    elem_nodes = spline_nodes[elem_node_ids, 1:end-1]
    return elem_nodes
end

function getCoefficientVectors( uspline::USpline )
    return uspline.coefficients.dense_coefficient_vectors
end

function getElementCoefficientVectorIds( uspline::USpline, elem_id )
    elem_idx = elemIdxFromElemId( uspline,  elem_id )
    return uspline.elements["element_blocks"][elem_idx]["coeff_vector_ids"] .+ 1
end

function getVertexConnectivity( uspline::USpline )
    return uspline.vertex_connectivity
end

function getElementExtractionOperator( uspline::USpline, elem_id )
    coeff_vectors = getCoefficientVectors( uspline )    
    coeff_vector_ids = getElementCoefficientVectorIds( uspline, elem_id )
    C = zeros( length( coeff_vector_ids ), length( coeff_vector_ids ) )
    for n = axes( coeff_vector_ids, 1 )
        C[n,:] = coeff_vectors[ coeff_vector_ids[n] ]
    end
    return C
end

function getGlobalExtractionOperator( uspline::USpline )
    num_elems = getNumElems( uspline )
    num_nodes = getNumNodes( uspline )
    num_bez_nodes = getNumBezierNodes( uspline )
    glob_extraction_operator = zeros( num_nodes, num_bez_nodes )
    for elem_idx = 1 : num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_node_ids = getElementNodeIds( uspline, elem_id )
        elem_bez_node_ids = getElementBezierNodeIds( uspline, elem_id )
        elem_extraction_operator = getElementExtractionOperator( uspline, elem_id )
        for i = axes( elem_bez_node_ids, 1 )
            I = elem_bez_node_ids[i]
            for j = axes( elem_bez_node_ids, 1 )
                J = elem_node_ids[j]
                glob_extraction_operator[J,I] = elem_extraction_operator[j, i]
            end
        end
    end
    return glob_extraction_operator
end

function getElementBezierNodes( uspline::USpline, elem_id )
    elem_nodes = getElementNodes( uspline, elem_id )
    C = getElementExtractionOperator( uspline, elem_id )
    element_bezier_node_coords = transpose( C ) * elem_nodes
    return element_bezier_node_coords
end

function getElementBezierVertices( uspline::USpline, elem_id )
    element_bezier_node_coords = getElementBezierNodes( uspline, elem_id )
    vertex_connectivity = getVertexConnectivity( uspline )
    vertex_coords = [ element_bezier_node_coords[1], element_bezier_node_coords[end] ]
    return vertex_coords
end

function getBezierNodes( uspline::USpline )
    bezier_nodes = []
    num_elems = getNumElems( uspline )
    for elem_idx = 1 : num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_bezier_nodes = getElementBezierNodes( uspline, elem_id )
        if isempty( bezier_nodes )
            bezier_nodes = elem_bezier_nodes
        else
            bezier_nodes = vcat( bezier_nodes, elem_bezier_nodes )
        end
    end
    bezier_nodes = uniquerows( bezier_nodes )
    return bezier_nodes
end

function getElementIdContainingPoint( uspline::USpline, point )
    num_elems = getNumElems( uspline )
    for elem_idx = 1 : num_elems
        elem_id = elemIdFromElemIdx( uspline, elem_idx )
        elem_domain = getElementDomain( uspline, elem_id )
        if ( ( point >= elem_domain[1] ) && ( point <= elem_domain[2] ) )
            return elem_id
        end
    end
end

function getNodeIdNearPoint( uspline::USpline, point )
    spline_nodes = getSplineNodes( uspline )[:,1]
    node_dist = sqrt( ( spline_nodes - point ) ^ 2.0 )
    return argminimum( node_dist )
end

function uniquerows( A::Matrix, atol = 0, rtol = 1e-12 )
    B = Matrix{ typeof( A[1] ) }( transpose( A[1,:] ) )
    for i = axes( A, 1 )
        is_row_eq = ones(Bool, size( B, 1 ) )
        for j = axes( B, 1 )
            is_row_eq[j] = isapprox( A[i,:], B[j,:], atol = atol, rtol = rtol )
        end
        if any( is_row_eq ) == false
            B = vcat( B, transpose( A[i,:] ) )
        end
    end
    return B
end