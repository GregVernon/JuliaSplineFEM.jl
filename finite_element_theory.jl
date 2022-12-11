### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 26bae7c1-3105-4734-9100-4008aed44ee0
import Pkg

# ╔═╡ 41990400-db45-44f7-800c-59b21759be30
Pkg.activate( "." )

# ╔═╡ ce426157-fcdb-4a03-a216-fddf1e2c1e87
using Gadfly

# ╔═╡ 0515b7ad-9911-4ed7-8b94-1fa5d441902b
using SpecialFunctions

# ╔═╡ c87841b6-dc3d-4f10-9f17-005001e9ead9
import JuliaSplineFEM

# ╔═╡ bd6e98dc-97b7-469e-8881-24fdc501ad44
md"""
### Target function
The function we want to approximate could be any arbitrary function, but here we limit ourselves to a simple univariate function
```math
f(x) = 
```
"""

# ╔═╡ d82dea7a-4783-4663-bae9-733808427ff7
f = (x) -> x * sin( 2*pi / 2 * x );

# ╔═╡ 620f7e8f-ed43-4ae5-8909-427e701229a5
md"""
## Build U-spline space
"""

# ╔═╡ 41c220fc-3408-43dc-ae34-9d88432bfb32
uspline = JuliaSplineFEM.readBEXT( "data/two_element_quadratic_bspline.json" );

# ╔═╡ f801cfdc-f41e-4945-b36a-9ccccd1428b6
d = JuliaSplineFEM.computeGalerkinApproximation( uspline, f );

# ╔═╡ 11840ea3-e6e1-4ac5-ad98-6d1e69b6f020
begin
	domain = JuliaSplineFEM.getDomain( uspline )
	x = LinRange( domain[1], domain[2], 1000 )
	plt = plot( )
	push!( plt, layer( f, domain[1], domain[2], color=[colorant"black"] ) )
	push!( plt, layer( (x)->JuliaSplineFEM.evaluateSolutionAt( uspline, d, x ), domain[1], domain[2] ) )
end

# ╔═╡ Cell order:
# ╠═26bae7c1-3105-4734-9100-4008aed44ee0
# ╠═41990400-db45-44f7-800c-59b21759be30
# ╠═ce426157-fcdb-4a03-a216-fddf1e2c1e87
# ╠═0515b7ad-9911-4ed7-8b94-1fa5d441902b
# ╠═c87841b6-dc3d-4f10-9f17-005001e9ead9
# ╠═bd6e98dc-97b7-469e-8881-24fdc501ad44
# ╠═d82dea7a-4783-4663-bae9-733808427ff7
# ╠═620f7e8f-ed43-4ae5-8909-427e701229a5
# ╠═41c220fc-3408-43dc-ae34-9d88432bfb32
# ╠═f801cfdc-f41e-4945-b36a-9ccccd1428b6
# ╠═11840ea3-e6e1-4ac5-ad98-6d1e69b6f020
