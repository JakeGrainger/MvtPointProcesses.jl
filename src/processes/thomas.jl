struct ThomasProcess{T1<:Real,T2<:Real,T3<:Real,G<:Geometry} <: PointProcess
    κ::T1
	σ::T2
    μ::T3
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

function sdf(p::ThomasProcess{T1, T2, T3, <:Geometry{2,S}},freq) where {T1,T2,T3,S}
    return p.κ*p.μ * (1+p.μ*exp(-p.σ^2* (2π)^2*norm(freq.coords)^2))
end