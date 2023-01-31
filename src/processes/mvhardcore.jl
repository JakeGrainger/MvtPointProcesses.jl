struct BivariateHardCore{D,T1,T2,T3,S,G<:Geometry{D,S}} <: PointProcess{D,2}
    rho_r::T1
    rho_c::T2
    r::T3
    geom::G
end

function rand(m::BivariateHardCore)
    X1 = rand(PoissonProcess(m.rho_r, m.geom))
    Y = rand(PoissonProcess(m.rho_c, m.geom))
    X2 = typeof(Y)[]
    for y in Y.items
        if !any(norm(x.coords .- y.coords) < m.r for x in X1.items)
            push!(X2, y)
        end
    end
    return (X1,PointSet(X2))
end