struct ThomasProcess{T,G<:Geometry} <: PointProcess
    κ::T
	σ::T
    μ::T
	geom::G
end

function rand(p::ThomasProcess)
    containing_box = boundingbox(p.geom)
    simulation_box = expandbox(containing_box,6p.σ)
    parents = rand(Poisson(p.κ, simulation_box))
    offspring = eltype(parents)[]
    for par ∈ parents
        n = rand(Poisson(p.μ))
        for _ ∈ 1:n
            push!(offspring, Point(rand.(Normal.(par.coords, p.σ))))
        end
    end
    return mask(PointSet(offspring), p.geom)
end