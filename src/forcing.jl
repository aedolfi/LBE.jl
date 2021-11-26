# want to include a force acting on the liquid. That is the place

"""
    function well(state::State, sys::SysConst)
creates velocity inflow from somewhere
"""
function well!(state:: State, sys:: SysConst)
    wellx = zeros(sys.Lx, sys.Ly)
    wellx[1,:] .+= 1e-6 + 1e-9 * (randn()-0.5)
    #wellx[sys.Lx,:] .-= 1e-6
    state.velx .+= wellx .+ 1e-9 * (randn()-0.5)
    state.vely .+= 1e-9 * (randn()-0.5)
    return nothing
end

function  sink!(state:: State, sys:: SysConst)
    state.velx[(sys.Lx-Int(round(3*sys.Lx/5))):(sys.Lx-Int(round(2*sys.Lx/5))), (sys.Ly-Int(round(3*sys.Ly/5))):(sys.Ly-Int(round(2*sys.Ly/5)))] .= 0.0
    state.vely[(sys.Lx-Int(round(3*sys.Lx/5))):(sys.Lx-Int(round(2*sys.Lx/5))), (sys.Ly-Int(round(3*sys.Ly/5))):(sys.Ly-Int(round(2*sys.Ly/5)))] .= 0.0
end