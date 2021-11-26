@testset "time time_loop" begin
    @testset "rho = 1 v = (0,0)" begin
        sys = SysConst{Float64}(Lx=50, Ly=50, Tmax=200)
        state = MyState(sys)
        result = time_loop(sys, state, verbose=false)
        mass = 0
        mass = sum(result.rho)
        @test all(mass ≈ 50^2)
    end
    @testset "rho=1, v=(1,0)" begin
        sys = SysConst{Float64}(Lx=50, Ly=50, Tmax=200)
        state = MyState(sys)
        initdis!(state, sys)
        result = time_loop(sys, state, verbose=false)
        @test all(mass ≈ 50^2)
end