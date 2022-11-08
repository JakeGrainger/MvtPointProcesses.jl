module MvtPointProcesses

using Distributions, Meshes, StaticArrays

include("processes/processes.jl")

export rand, PoissonProcess, CoxProcess, coxprocess

end