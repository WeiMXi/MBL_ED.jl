using MBL_ED
using Test
using JLD

@testset "MBL_ED.jl" begin
    tasks = give_tasks( tol_n = 10
    , up = 5
    , J = 1.0
    , Jp = 1.0
    , Jz = 1.0
    , W = 0.5:0.5:8.0
    , k = (sqrt(5) - 1)/2
    , samples = 10
    )
    t_data = load("./test_data.jld")
    r1_t = cal_one_task(t_data["task"])
    r2_t = primitive_run(tol_n = 10, up = div(10, 2), W = 1.0:0.5:2.0, samples = 10)
    @test tasks[2].H == t_data["task"].H
    @test abs(r1_t[:, 4][1] - t_data["result"][4]) < 0.01
    @test abs(r1_t[:, 5][1] - t_data["result"][5]) < 0.01
    @test 0 < sum(r2_t.W)/length(r2_t.W) < Inf
    @test 0 < sum(r2_t.r)/length(r2_t.r) < Inf
    # @test abs.( cal_one_task(t_data["task"])[4:5] .- t_data["result"][4:5] ) .< 0.01
end
