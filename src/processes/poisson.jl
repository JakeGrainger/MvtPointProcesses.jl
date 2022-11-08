struct IntensityGrid{D,T1,T2}
	ρ::Array{T1,D}
	grid::CartesianGrid{D,T2}
	function IntensityGrid(ρ::Array{T1,D}, grid::CartesianGrid{D,T2}) where {D,T1,T2}
		@assert size(ρ) == size(grid)
		new{D,T1,T2}(ρ, grid)
	end
end
struct Intensity{T<:Real,F<:Union{Function,IntensityGrid}}
	ρ::F
	ρ₀::T
end

struct PoissonProcess{T<:Union{Real,Intensity},G<:Geometry} <: PointProcess
	ρ::T
	geom::G
end

function (ig::IntensityGrid)(ξ)
	ξ ∈ g || error("point is not in domain of intensity grid.")
	return findfirst(ξ in g for g in ig.grid)
end
Intensity(ρ::IntensityGrid) = Intensity(ρ, maximum(ρ))
Intensity(ρ::Function, mesh::Mesh) = maximum(ρ(centroid(m)) for m in mesh)

Base.maximum(g::IntensityGrid) = Base.maximum(g.ρ)

function rand(p::PoissonProcess{<:Real,<:Geometry})
	grid = boundinggrid(p.geom)
	N = rand(Poisson(p.ρ * measure(grid)))
	X = PointSet(collect(sample(grid, HomogeneousSampling(N))))
	return mask(X,p.geom)
end

function rand(p::PoissonProcess{<:Intensity,<:Geometry})
	X = Base.rand(PoissonProcess(p.ρ.ρ₀, p.geom))
	thin(X, ξ->p.ρ.ρ(ξ)/p.ρ₀)
end