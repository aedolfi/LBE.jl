@testset "Bounce back Boundaries" begin
    @testset "boundarie clasification" begin
        sys = LBE.SysConstWithBound{Float64}(Lx=50, Ly=50, Tmax=200)
        obs = zeros(sys.Lx, sys.Ly)
        obsleft = zeros(sys.Lx, sys.Ly)
        obsright = zeros(sys.Lx, sys.Ly)
        obsup = zeros(sys.Lx, sys.Ly)
        obsdown = zeros(sys.Lx, sys.Ly)
        view(obs, :, 1) .= 1
        view(obs, 1,: ) .= 1
        view(obs, :, sys.Ly) .= 1
        view(obs, sys.Lx,: ) .= 1
        view(obsleft, :, 1) .= 1
        view(obsup, 1,: ) .= 1
        view(obsright, :, sys.Ly) .= 1
        view(obsdown, sys.Lx,: ) .= 1
        sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
        @test all(sys.obsdown  .= obsdown)
        @test all(sys.obsleft .= obsleft)
        @test all(sys.obsright  = obsright)
        @test all(sys.obsup .= obsup)
    end
    @testset "square" begin
        sys = LBE.SysConstWithBound{Float64}(Lx=50, Ly=50, Tmax=200)
        obs = zeros(sys.Lx, sys.Ly)
        view(obs, :, 1) .= 1
        view(obs, 1,: ) .= 1
        view(obs, :, sys.Ly) .= 1
        view(obs, sys.Lx,: ) .= 1
        sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
        state = LBE.MyState(sys)
        result = LBE.time_loop_channel(sys, state, verbose=false)
        mass = 0
        mass = sum(result.rho)
        @test all(mass ≈ 50^2)
    end
    @testset "upright channel" begin
        sys = LBE.SysConstWithBound{Float64}(Lx=50, Ly=50, Tmax=200)
        obs = zeros(sys.Lx, sys.Ly)
        view(obs, 1,: ) .= 1
        view(obs, sys.Lx,: ) .= 1
        sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
        state = LBE.MyState(sys)
        result = LBE.time_loop_channel(sys, state, verbose=false)
        mass = 0
        mass = sum(result.rho)
        @test all(mass ≈ 50^2)
    end
    @testset "channel" begin
        sys = LBE.SysConstWithBound{Float64}(Lx=50, Ly=50, Tmax=200)
        obs = zeros(sys.Lx, sys.Ly)
        view(obs, :, 1) .= 1
        view(obs, :, sys.Ly) .= 1
        sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
        state = LBE.MyState(sys)
        result = LBE.time_loop_channel(sys, state, verbose=false)
        mass = 0
        mass = sum(result.rho)
        @test all(mass ≈ 50^2)
    end
    @testset "channel rho = 1, v =(1,0)" begin
        sys = LBE.SysConstWithBound{Float64}(Lx=50, Ly=50, Tmax=200)
        obs = zeros(sys.Lx, sys.Ly)
        view(obs, 1, :) .= 1
        view(obs, sys.Lx, : ) .= 1
        sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
        state = LBE.MyStateX(sys)
        LBE.initdis!(state, sys)
        result = LBE.time_loop_channel(sys, state, verbose = false)
        @test all(mass ≈ 50^2)
end