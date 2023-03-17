@testset "sdf tests" begin
    model = PoissonProcess(0.01, Box(Point(0,0), Point(1,1)))
    @test sdf(model, 0.1) == 0.01
end