module MvtPointProcesses

using Distributions, Meshes, StaticArrays
import Base: rand

include("utils.jl")
include("processes/processes.jl")

export PoissonProcess, CoxProcess, coxprocess

end