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
    B = Bext( is_rational, num_nodes, num_vertices, num_elems, dim, num_coeff_vecs, patch_id, elements, block_elem_map, vertex_connectivity, nodes, coefficients)
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

struct Bext
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

