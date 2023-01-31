struct ShiftedProcess{D,P,M<:PointProcess{D,P},T<:NTuple{P,NTuple{D,Real}}} <: PointProcess{D,P}
    original::M
    shift::T
end
shift(p::PointProcess, s) = ShiftedProcess(p, s)
shift(p::PointSet{D,<:Real}, s::NTuple{D,Real}) where {D} = PointSet([Point(x.coords .+ s) for x in p.items])

rand(p::ShiftedProcess{D,1}) where {D} = rand(p)
function rand(p::ShiftedProcess{D,P}) where {D,P}
    X = rand(expand_geom(p.original, max_shift(p.shift)))
    Y = shift.(X,p.shift)
    return mask.(Y, Ref(p.original.geom))
end

max_shift(s) = maximum(x->maximum(abs,x), s)