using LBE, Test

@testset "LBE.jl" begin
    println("boundaries")
    include("boundaries'jl")
    println("time_loop")
    include("time_loop.jl")
end