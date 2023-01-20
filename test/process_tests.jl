function test_point_process(model::MvtPointProcesses.PointProcess{D,P}) where {D,P}
    @testset "Tests for model: $model." begin
        X = rand(model)
        @test X isa expected_sim_output(model)
    end
end
expected_sim_output(::MvtPointProcesses.PointProcess{D,P}) where {D,P} = NTuple{P,<:PointSet}
expected_sim_output(::MvtPointProcesses.PointProcess{D,1}) where {D} = PointSet
function test_box()
    return Box(Point(0,0),Point(100,100))
end

@testset "Tests for Poisson" begin
    test_point_process(PoissonProcess(0.01, test_box()))
    test_point_process(PoissonProcess(MvtPointProcesses.Intensity(x->0.01, 0.02), test_box()))
end

@testset "Tests for Thomas" begin
    test_point_process(ThomasProcess(0.01, 4, 0.05, test_box()))
    test_point_process(ThomasProcess(0.01, (4,3), (0.05,0.06), test_box()))
end

@testset "Tests for MaternI" begin
    test_point_process(MaternIProcess(0.01, 1, test_box()))
end

@testset "Tests for MaternII" begin
    test_point_process(MaternIIProcess(0.01, 1, test_box()))
end