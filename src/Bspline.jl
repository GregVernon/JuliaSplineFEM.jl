struct SplineSpace
    degree::Int
    vertices::Array
    continuity::Array
end

function evalBspline( space::SplineSpace, variate )
    degree = space.degree
    kV = knotVectorFromSplineSpace( space )
    N = []
    for p = 0:degree
        push!( N, zeros( length( kV ) - ( p + 1 ) ) )
        for ii = 1 : ( length( kV ) - ( p + 1 ) )
            if p == 0
                if kV[ii] <= variate && variate < kV[ii + 1]
                    N[p+1][ii] = 1
                else
                    N[p+1][ii] = 0
                end
            else
                divisor_1 =  kV[ii + p] - kV[ii]
                if divisor_1 == 0
                    term_1 = 0
                else
                    term_1 = ( variate - kV[ii] ) / ( kV[ii + p] - kV[ii] )
                end

                divisor_2 = kV[ii + p + 1] - kV[ii + 1]
                if divisor_2 == 0
                    term_2 = 0
                else
                    term_2 = ( kV[ii + p + 1] - variate ) / ( kV[ii + p + 1] - kV[ii + 1]);
                end
                N[p + 1][ii] = term_1 * N[p][ii] + term_2 * N[p][ii + 1];
            end
        end
    end
    return N[end]
end

function knotVectorFromSplineSpace( space::SplineSpace )
    knotVector = []
    for ii = 1 : length(space.vertices)
        append!( knotVector, repeat( [ space.vertices[ii] ], space.degree - space.continuity[ii] ) );
    end
    knotVector = vcat( kv... )
    return knotVector
end

function globalBSplineBasisIdFromLocalElemId( space::SplineSpace, elem_id, local_fun_id )
    global_fun_id = 0;
    for e = 1 : elem_id
        elem_left_continuity = splineSpace.continuity[ e ];
        if e == 1 && elem_id == 1
            global_fun_id = local_fun_id;
        elseif e == 1
            global_fun_id = global_fun_id + splineSpace.degree;
        elseif e == elem_id
            global_fun_id = global_fun_id + ( local_fun_id - elem_left_continuity );
        else
            global_fun_id = global_fun_id + ( splineSpace.degree - elem_left_continuity );
        end
    end
    return global_fun_id
end
    
function elementGlobalBSplineBasisIdFromLocalElemId( space::SplineSpace, elem_id )
degree = splineSpace.degree;
elem_global_fun_ids = zeros( degree + 1 );
for ii = 1 : length( elem_global_fun_ids )
    elem_global_fun_ids[ii] = globalSplineBasisIdFromLocalElemId( splineSpace, elem_id, ii );
end
return elem_global_fun_ids
end