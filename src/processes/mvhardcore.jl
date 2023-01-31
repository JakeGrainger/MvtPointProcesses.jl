struct BivariateHardCoreProcess{D,T1,T2,T3,S,G<:Geometry{D,S}} <: PointProcess{D,2}
    rho_r::T1
    rho_c::T2
    r::T3
    geom::G
end

function rand(m::BivariateHardCoreProcess)
    bigbox = expandbox(boundingbox(m.geom), m.r)
    X1 = rand(PoissonProcess(m.rho_r, bigbox))
    Y = rand(PoissonProcess(m.rho_c, bigbox))
    X2 = eltype(Y)[]
    for y in Y.items
        if !any(norm(x.coords .- y.coords) < m.r for x in X1.items)
            push!(X2, y)
        end
    end
    return filter.((X1,PointSet(X2)), Ref(m.geom))
end