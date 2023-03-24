struct IntensityGrid{D,T1,T2}
	λ::Array{T1,D}
	grid::CartesianGrid{D,T2}
	function IntensityGrid(λ::Array{T1,D}, grid::CartesianGrid{D,T2}) where {D,T1,T2}
		@assert size(λ) == size(grid)
		new{D,T1,T2}(λ, grid)
	end
end
struct Intensity{T<:Real,F<:Union{Function,IntensityGrid}}
	λ::F
	λ₀::T
end

"""
	PoissonProcess(λ, geom)

A Poisson process with intensity λ to be sampled on geometry `geom`.
The intensity can be `Real`, in which case the process is stationary, or `Intensity`.
"""
struct PoissonProcess{D,T,S<:Union{Real,Intensity},G<:Geometry{D,T}} <: PointProcess{D,1}
	λ::S
	geom::G
end

function (g::IntensityGrid)(ξ)
	gmin = minimum(g.grid).coords
    gΔ = g.grid.spacing
    gn = g.grid.topology.dims
	all(gmin.≤ ξ.coords .≤gmin.+gΔ.*gn) || error("point is not in domain of intensity grid.")
    ind = CartesianIndex(Tuple(min.(floor.(Int, (ξ.coords.-gmin)./gΔ).+1, gn) ))
    return g.λ[ind]
end
Intensity(λ::IntensityGrid) = Intensity(λ, maximum(λ))
Intensity(λ::Function, mesh::Mesh) = maximum(λ(centroid(m)) for m in mesh)

Base.maximum(g::IntensityGrid) = Base.maximum(g.λ)

function rand(p::PoissonProcess{D,T,<:Real,<:Geometry{D,T}}) where {D,T}
	grid = boundingbox(p.geom)
	N = rand(Poisson(p.λ * measure(grid)))
	U = ntuple(d-> Uniform(minimum(grid).coords[d], maximum(grid).coords[d]), Val{D}())
	X = PointSet([Point(rand.(U)) for _ in 1:N])
	return mask(X,p.geom)
end

function rand(p::PoissonProcess{D,T,<:Intensity,<:Geometry}) where {D,T}
	X = Base.rand(PoissonProcess(p.λ.λ₀, p.geom))
	thin(X, ξ->p.λ.λ(ξ)/p.λ.λ₀)
end

sdf(p::PoissonProcess{D,T,<:Real,<:Geometry},freq) where {D,T} = p.λ
sdf(p::PoissonProcess{D,T,<:Intensity,Geometry},freq) where {D,T} = error("Process is not stationary, so sdf not defined.")