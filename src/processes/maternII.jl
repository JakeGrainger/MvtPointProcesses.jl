struct MaternIIProcess{T1<:Real,T2<:Real,G<:Geometry} <: PointProcess
	κ::T1
    r::T2
	geom::G
end

function rand(p::MaternIIProcess)
    proposals = rand(PoissonProcess(p.κ, expandbox(boundingbox(p.geom), p.r))) # generate on a bounding box inflated by the inhibition radius
    marks = rand(length(proposals))
    X = eltype(proposals.items)[]
    for (m1,p) in zip(marks,proposals.items) #  note this could be made more efficient
        if sum(m1 ≤ m2 && norm(p.coords.-q.coords) < r for (m2,q) in zip(marks,proposals.items)) == 1 # one of these is distance to self
            push!(X,p)
        end
    end
    return mask(PointSet(X),p.geom)
end