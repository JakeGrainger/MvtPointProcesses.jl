function boundinggrid(geom::Geometry{D,T},grid_res::NTuple{D,Int}=ntuple(d->1000, Val{D}())) where {D,T}
    box = boundingbox(geom)
    return CartesianGrid(box.min,box.max,dims=grid_res)
end

function expandbox(box::Box, pad::Number)
	return Box(Point(box.min.coords .- pad), Point(box.max.coords .+ pad))
end

function mask(X::PointSet,geom::Geometry)
    return PointSet([x for x ∈ X if x ∈ geom])
end

format_pp_output(X::NTuple{N,T}) where {N,T<:PointSet} = X
format_pp_output(X::NTuple{1,T}) where {T<:PointSet} = X[1]

function thin(X::PointSet, retention_probability)
	thinned = eltype(X)[]
	for ξ ∈ X
		if rand() ≤ retention_probability(ξ)
			push!(thinned, ξ)
		end
	end
	return PointSet(thinned)
end

function expand_geom(p::P,r) where {P<:PointProcess}
	args = map(name->name==:geom ? expandbox(boundingbox(p.geom), r) : getfield(p, name), fieldnames(P))
	P(args...)
end