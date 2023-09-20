"""
    MaternIProcess(κ, r, geom)
A Matern hard core process of type I.
Proposal points from a poisson process are first generated with intensity `κ`.
Each point is included in the final process if no other point lies within a distance of `r` of the point in question.

# Parameters
- `κ`: The intensity of the proposal locations.
- `r`: The repulsion radius.
- `geom`: The geometry on which we wish to sample the process.
"""    
struct MaternIProcess{D,T,T1<:Real,T2<:Real,G<:Geometry{D,T}} <: PointProcess{D,1}
	κ::T1
    r::T2
	geom::G
end

function rand(model::MaternIProcess)
    proposals = rand(PoissonProcess(model.κ, expandbox(boundingbox(model.geom), model.r))) # generate on a bounding box inflated by the inhibition radius
    X = eltype(parent(proposals))[]
    for p in parent(proposals) #  note this could be made more efficient
        if sum(norm(p.coords.-q.coords) < model.r for q in parent(proposals)) == 1 # one of these is distance to self
            push!(X,p)
        end
    end
    return mask(PointSet(X), model.geom)
end