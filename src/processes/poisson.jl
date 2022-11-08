abstract type PoissonProcess <: PointProcess end

struct HomogeneousPoissonProcess{T<:Real} <: PoissonProcess
	ρ::T
end

struct InHomogeneousPoissonProcess{T<:Real,F<:Function} <: PoissonProcess
	ρ₀::T
	ρ::F
end

function simulate(p::HomogeneousPoissonProcess,B::Box)
	N = rand(Poisson(p.ρ * prod(B.a)))
	dist = Uniform.(0,B.a)
	points = [rand.(dist) for i in 1:N]
	return points
end

function simulate(p::InHomogeneousPoissonProcess,B::Box)
	points = simulate(HomogeneousPoissonProcess(p.ρ₀),B)
	thinned = eltype(points)[]
	for ξ ∈ points
		if rand() ≤ p.ρ(ξ)/p.ρ₀
			push!(thinned, ξ)
		end
	end
	return thinned
end