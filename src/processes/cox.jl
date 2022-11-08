"""
    CoxProcess(Λ,link, geom)

Cox process with random intensity field `Λ` and link function `link` over geometry `geom`.

The intensity field should have a method rand which allows for the generation of a random field.
"""
struct CoxProcess{D,F<:Function,G} <: PointProcess
    Λ::D
    link::F
    geom::G
end

"""
    CoxProcess(geom::Geometry, fieldtype, link; grid_res)

Produces a CoxProcess over the specified geometry `geom` with a field of type `fieldtype` and link function `link`.

Note that `fieldtype` should be a function which produces a random field when provided with a mesh.

Optionally specify `grid_res` to get a better quality simulation.
"""
function CoxProcess(fieldtype, link, geom::Geometry{D,T}; grid_res::NTuple{D,Int}=ntuple(d->1000,val{D}())) where {D}
    box = boundingbox(geom)
    grid = CartesianGrid(box.min,box.max,dims=grid_res)
    CoxProcess(fieldtype(grid), link, geom)
end
CoxProcess(fieldtype, link, geom::Geometry{D,T}; grid_res::Int) = CoxProcess(fieldtype, link, geom::Geometry{D,T}; grid_res = ntuple(d->grid_res,Val{D}())) where {D}

function rand(c::CoxProcess)
    intensity = c.link(rand(c.Λ)) # generate intensity field
    X_grid = rand.(PoissonProcess.(indensity)) # generate inhomogeneous Poisson processes
    X = PointSet([x for x in X_grid if x ∈ c.geom]) # keep only points in the domain
    return (X=X,intensity=intensity)
end