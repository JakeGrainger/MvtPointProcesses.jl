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
    grid = boundinggrid(geom)
    CoxProcess(fieldtype(grid), link, geom)
end
CoxProcess(fieldtype, link, geom::Geometry{D,T}; grid_res::Int) = CoxProcess(fieldtype, link, geom::Geometry{D,T}; grid_res = ntuple(d->grid_res,Val{D}())) where {D}

function Base.rand(c::CoxProcess)
    intensity = rand(c.Λ) # generate intensity field
    return _cox_rand(intensity, c)
end

_cox_rand(intensity::Array{SVector{P,T},D}, c::CoxProcess) where {P,D,T} = _cox_rand(ntuple(p->getindex.(intensity,p), Val{P}()), c)
_cox_rand(intensity::NTuple{P,Array{T,D}}, c::CoxProcess) where {P,D,T} = _cox_rand.(intensity, Ref(c))
function _cox_rand(intensity::Array{T,D}, c::CoxProcess) where {D,T}
    transformed_intensity = c.link.(intensity)
    X = rand(PoissonProcess(transformed_intensity, c.geom)) # generate inhomogeneous Poisson processes
    return (X=X,intensity=transformed_intensity)
end
