struct MultivariateHardCoreProcess{D,P,T1,T2,T3,S,G<:Geometry{D,S}} <: PointProcess{D,P}
    rho_r::T1
    rho_c::T2
    r::T3
    geom::G
    MultivariateHardCoreProcess(rho_r::T1, rho_c::T2, r::T3, geom::G, ::Val{P}) where {T1,T2,T3,G<:Geometry{D,S},P} where {D,S} = new{D,P,T1,T2,T3,S,G}(rho_r, rho_c, r, geom)
    MultivariateHardCoreProcess(rho_r, rho_c::Real, r, geom) = MultivariateHardCoreProcess(rho_r, (rho_c,), r, geom, Val{2}())
    MultivariateHardCoreProcess(rho_r, rho_c::NTuple{Q,Real}, r, geom) where {Q} = MultivariateHardCoreProcess(rho_r, rho_c, r, geom, Val{Q+1}())
end

const BivariateHardCoreProcess{D,T1,T2,T3,S,G} = MultivariateHardCoreProcess{D,2,T1,T2,T3,S,G}
BivariateHardCoreProcess(args...) = MultivariateHardCoreProcess(args...)

function rand(m::MultivariateHardCoreProcess)
    bigbox = expandbox(boundingbox(m.geom), m.r)
    X1 = rand(PoissonProcess(m.rho_r, bigbox))
    Y = rand.(PoissonProcess.(m.rho_c, Ref(m.geom)))
    X2 = inhibit.(Y, Ref(X1), m.r)
    return (mask(X1, m.geom), X2...)
end


"""
    inhibit(Y::PointSet,X::PointSet,r::Real)

Inhibits Y based on X with radius r.
"""
function inhibit(Y::PointSet,X::PointSet,r::Real)
    Y_inhib = eltype(Y)[]
    for y in Y.items
        if !any(norm(x.coords .- y.coords) < r for x in X.items)
            push!(Y_inhib, y)
        end
    end
    return PointSet(Y_inhib)
end