module JuliaSplineFEM

import SpecialMatrices
import LinearAlgebra
import Polynomials
import ForwardDiff
import JSON
using Memoize

include( "Bext.jl" )
include( "Basis.jl" )
include( "Quadrature.jl" )
include( "Galerkin.jl" )

end
