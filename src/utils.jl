function boundinggrid(geom::Geometry{D,T},grid_res::NTuple{D,Int}=ntuple(d->1000, Val{D}())) where {D,T}
    box = boundingbox(geom)
    return CartesianGrid(box.min,box.max,dims=grid_res)
end

function mask(X::PointSet,geom::Geometry)
    return PointSet([x for x ∈ X if x ∈ geom])
end

function thin(X, retention_probability)
	thinned = eltype(X)[]
	for ξ ∈ X
		if rand() ≤ retention_probability(ξ)
			push!(thinned, ξ)
		end
	end
	return thinned
end