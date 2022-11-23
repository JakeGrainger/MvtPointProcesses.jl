module MvtPointProcesses

using Distributions, Meshes, StaticArrays, CirculantEmbedding
import Base: rand

include("utils.jl")
include("processes/processes.jl")

export PoissonProcess, CoxProcess, coxprocess

end