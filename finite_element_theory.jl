import Pkg
Pkg.activate( "." )
using Gadfly
using SpecialFunctions
import JuliaSplineFEM

uspline = JuliaSplineFEM.readBEXT( "data/two_element_quadratic_bspline.json" )
f = (x) -> sin( pi / 2 * x )
d = JuliaSplineFEM.computeGalerkinApproximation( uspline, f )

domain = JuliaSplineFEM.getDomain( uspline )
x = LinRange( domain[1], domain[2], 1000 )
plt = plot( )
push!( plt, layer( f, domain[1], domain[2], color=[colorant"black"] ) )
push!( plt, layer( (x)->JuliaSplineFEM.evaluateSolutionAt( uspline, d, x ), domain[1], domain[2] ) )
img = SVG("function.svg", 6inch, 4inch)
draw(img, plt)