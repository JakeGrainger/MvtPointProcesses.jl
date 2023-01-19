struct MaternIProcess{T1<:Real,T2<:Real,G<:Geometry} <: PointProcess
	κ::T1
    r::T2
	geom::G
end

function rand(p::MaternIProcess)
    proposals = rand(PoissonProcess(p.κ, expandbox(boundingbox(p.geom), p.r))) # generate on a bounding box inflated by the inhibition radius
    X = eltype(proposals.items)[]
    for p in proposals.items #  note this could be made more efficient
        if sum(norm(p.coords.-q.coords) < r for q in proposals.items) == 1 # one of these is distance to self
            push!(X,p)
        end
    end
    return mask(PointSet(X),p.geom)
end