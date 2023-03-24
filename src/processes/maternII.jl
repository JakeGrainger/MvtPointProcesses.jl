"""
    MaternIIProcess

A Matern hard core process of type II.
Generates a marked proposal process with intensity `κ` and marks with iid Unif(0,1).
Points are included if they have the largest mark of all points within a radius `r` of the point in question.

# Parameters
- `κ`: The intensity of the proposal locations.
- `r`: The repulsion radius.
- `geom`: The geometry on which we wish to sample the process.
"""
struct MaternIIProcess{D,T,T1<:Real,T2<:Real,G<:Geometry{D,T}} <: PointProcess{D,1}
	κ::T1
    r::T2
	geom::G
end

function rand(model::MaternIIProcess)
    proposals = rand(PoissonProcess(model.κ, expandbox(boundingbox(model.geom), model.r))) # generate on a bounding box inflated by the inhibition radius
    marks = rand(length(proposals))
    X = eltype(proposals.items)[]
    for (m1,p) in zip(marks,proposals.items) #  note this could be made more efficient
        if sum(m1 ≤ m2 && norm(p.coords.-q.coords) < model.r for (m2,q) in zip(marks,proposals.items)) == 1 # one of these is distance to self
            push!(X,p)
        end
    end
    return mask(PointSet(X), model.geom)
end