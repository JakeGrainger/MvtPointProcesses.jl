module MvtPointProcesses

using Distributions, Meshes, StaticArrays

include("utils.jl")
include("processes/processes.jl")

export rand, PoissonProcess, CoxProcess, coxprocess

end