struct ThomasProcess{D,P,T,T1<:Real,T2<:NTuple{P,<:Real},T3<:NTuple{P,<:Real},G<:Geometry{D,T}} <: PointProcess{D,P}
    κ::T1
	σ::T2
    μ::T3
	geom::G
end
ThomasProcess(κ,σ::Real,μ::Real,geom) = ThomasProcess(κ,(σ,),(μ,),geom)

function rand(p::ThomasProcess)
    containing_box = boundingbox(p.geom)
    simulation_box = expandbox(containing_box,6maximum(p.σ))
    parents = rand(PoissonProcess(p.κ, simulation_box))
    offspring = thomas_offspring.(p.σ, p.μ, Ref(parents))
    return format_pp_output(mask.(offspring, Ref(p.geom)))
end

function thomas_offspring(σ::Real, μ::Real, parents)
    offspring = eltype(parents)[]
    for par ∈ parents
        n = rand(Poisson(μ))
        for _ ∈ 1:n
            push!(offspring, Point(rand.(Normal.(par.coords, σ))))
        end
    end
    return PointSet(offspring)
end

function sdf(p::ThomasProcess{2,1,T,T1,T2,T3,G},freq) where {T,T1,T2,T3,G}
    return p.κ*p.μ * (1+p.μ*exp(-p.σ^2* (2π)^2*norm(freq.coords)^2))
end