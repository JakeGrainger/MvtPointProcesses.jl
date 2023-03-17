using MvtPointProcesses
import MvtPointProcesses: SVector, SMatrix
using Test

@testset "MvtPointProcesses.jl" begin
    include("process_tests.jl")
    include("sdf_test.jl")
end
