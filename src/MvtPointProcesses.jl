module MvtPointProcesses

using Distributions, Meshes, StaticArrays

include("utils.jl")
include("processes/processes.jl")

export simulate, PoissonProcess, CoxProcess, coxprocess

end