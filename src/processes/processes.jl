abstract type PointProcess{D,P} end

include("poisson.jl")
include("cox.jl")
include("thomas.jl")
include("maternI.jl")
include("maternII.jl")