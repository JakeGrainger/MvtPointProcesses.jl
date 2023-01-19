function test_point_process(model::MvtPointProcesses.PointProcess)
    @testset "Tests for model: $model." begin
        X = rand(model)
        @test X isa PointSet
    end
end
function test_box()
    return Box(Point(0,0),Point(100,100))
end

@testset "Tests for Poisson" begin
    test_point_process(PoissonProcess(0.01, test_box()))
    test_point_process(PoissonProcess(Intensity(x->0.01, 0.02), test_box()))
end

@testset "Tests for Thomas" begin
    test_point_process(ThomasProcess(0.01, 4, 0.05, test_box()))
end

@testset "Tests for MaternI" begin
    test_point_process(MaternIProcess(0.01, 1, test_box()))
end

@testset "Tests for MaternII" begin
    test_point_process(MaternIIProcess(0.01, 1, test_box()))
end