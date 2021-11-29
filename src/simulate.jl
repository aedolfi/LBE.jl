"""
    function time_loop(sys:: SysConst, state::State: verbose=true)
Runs the simulation
"""
function time_loop(sys:: SysConst, state::State; verbose=true)
    for t in 1:sys.Tmax
        if t % sys.tdump==0
            mass = 0.0
            mass = sum(state.rho)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        #well!(state, sys)
        #sink!(state,sys)
        LBE.equilibrium!(state, sys)
        LBE.BGKandStream!(state, sys)
        LBE.moments!(state, sys)
    end
    return state
end


"""
    function time_loop_channel(sys, state, verbose)
Runs with channel boundary where every right pointing velocity in obsright
 and left pointing velocity in obsleft is forgotten and any up/down pointing velocity in obsup/obsdown 
 is reflected
"""
function time_loop_channel(sys:: SysConstWithBound, state::StateWithBound; verbose=true)
    for t in 1:sys.Tmax
        if t % sys.tdump==0
            mass = 0.0
            mass = sum(state.rho)
            if verbose
                println("Time step $t mass is $(round(mass, digits=3))")
            end
        end
        #well!(state, sys)
        #sink!(state,sys)
        LBE.equilibrium!(state, sys)
        LBE.BGKandStreamChannel!(state, sys)
        LBE.moments!(state, sys)
    end
    return state
end






function main()
    sys = LBE.SysConst{Float64}(Tmax=1000)
    state = LBE.MyStateX(sys)
    LBE.initdis!(state, sys)
    result = LBE.time_loop(sys, state)
    plot(heatmap(result.velsq),heatmap(result.rho))
end


function channel()
    sys = LBE.SysConstWithBound{Float64}()
    obs = zeros(sys.Lx, sys.Ly)
    #view(obs, 20:21 , : ).=1
    for i in 2:8, j in 2:3
        obs[i,j] = 1
    end
    sys = LBE.SysConstWithBound{Float64}(Tmax=1000, obs=obs)
    state = LBE.MyState(sys)
    #LBE.initdis!(state, sys)
    result = LBE.time_loop_channel(sys, state)
    plot(heatmap(result.velsq),heatmap(result.rho))
end