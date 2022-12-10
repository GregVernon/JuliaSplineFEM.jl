module JuliaSplineFEM

import SpecialMatrices
import LinearAlgebra
import Polynomials
import ForwardDiff
import JSON
using Memoize

include( "Basis.jl" )
include( "Quadrature.jl" )
include( "Bext.jl" )

end
