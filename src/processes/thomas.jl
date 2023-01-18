struct ThomasProcess{T,G<:Geometry} <: PointProcess
    κ::T
	σ::T
    μ::T
	geom::G
end

function rand(p::ThomasProcess)
    containing_box = boundingbox(p.geom)
    simulation_box = expandbox(containing_box,6p.σ)
    parents = rand(PoissonProcess(p.κ, simulation_box))
    offspring = eltype(parents)[]
    for par ∈ parents
        n = rand(Poisson(p.μ))
        for _ ∈ 1:n
            push!(offspring, Point(rand.(Normal.(par.coords, p.σ))))
        end
    end
    return mask(PointSet(offspring), p.geom)
end

function sdf(p::ThomasProcess{T1, <:Geometry{2,T2}},k) where {T1,T2}
    return p.λ*p.μ * (1+p.μ*exp(-p.σ^2* (2π)^2*norm(k.coords)^2))
end