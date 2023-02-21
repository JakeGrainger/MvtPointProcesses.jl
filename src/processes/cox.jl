"""
    CoxProcess(Λ,link, geom)

Cox process with random intensity field `Λ` and link function `link` over geometry `geom`.

The intensity field should have a method rand which allows for the generation of a random field.
"""
struct CoxProcess{D,P,T,S<:RandomField{D,P},F<:Function,G<:Geometry{D,T}} <: PointProcess{D,P}
    Λ::S
    link::F
    geom::G
end

"""
    CoxProcess(geom::Geometry, fieldtype, link; grid_res)

Produces a CoxProcess over the specified geometry `geom` with a field of type `fieldtype` and link function `link`.

Note that `fieldtype` should be a function which produces a random field when provided with a mesh.

Optionally specify `grid_res` to get a better quality simulation.
"""
function coxprocess(fieldtype, link, geom::Geometry{D,T}, grid_res::NTuple{D,Int}) where {D,T}
    grid = boundinggrid(geom,grid_res)
    CoxProcess(fieldtype(grid), link, geom)
end
coxprocess(fieldtype, link, geom::Geometry{D,T}, grid_res::Int=1000) where {D,T} = coxprocess(fieldtype, link, geom, ntuple(d->grid_res,Val{D}()))

function Base.rand(c::CoxProcess)
    intensity = Base.rand(c.Λ) # generate intensity field
    return process_fields(_cox_rand(striplatent(intensity), c), getlatent(intensity))
end

striplatent(x) = x
striplatent(x::NamedTuple{(:rf,:latent), T}) where {T} = x.rf
getlatent(x) = nothing
getlatent(x::NamedTuple{(:rf,:latent), T}) where {T} = x.latent

process_fields(x, ::Nothing) = x
process_fields(x, y) = (x..., latent=y)

_cox_rand(intensity::Array{SVector{P,T},D}, c::CoxProcess) where {P,D,T} = _cox_rand(ntuple(p->getindex.(intensity,p), Val{P}()), c)
_cox_rand(intensity::NTuple{P,Array{T,D}}, c::CoxProcess) where {P,D,T} = _cox_rand.(intensity, Ref(c))
function _cox_rand(intensity::Array{T,D}, c::CoxProcess) where {D,T}
    transformed_intensity = Intensity(IntensityGrid(c.link.(intensity), getmesh(c.Λ)))
    X = Base.rand(PoissonProcess(transformed_intensity, c.geom)) # generate inhomogeneous Poisson processes
    return (X=X,intensity=transformed_intensity.ρ.ρ)
end

function Distributions.mean(c::CoxProcess{D,P,T,S,typeof(exp),G}) where {D,P,T,S,G}
    return exp(mean(c.Λ) + var(c.Λ)/2)
end

function Distributions.cor(c::CoxProcess{D,P,T,S,typeof(exp),G}, h) where {D,P,T,S,G}   
    return exp(cov(c.Λ, h))
end

function Distributions.cov(c::CoxProcess{D,P,T,S,typeof(exp),G}, h) where {D,P,T,S,G}
    return mean(c)^2 * cor(c, h)
end

function approximate_cov(c::CoxProcess{D,P,T,S,typeof(exp),G}, lags) where {D,P,T,S,G}
    cov_field = approx_cov(c.Λ, lags)
    pp_mean = exp(mean(c.Λ) + cov_field[1]/2)
    return pp_mean^2 .* exp.(cov_field), pp_mean
end