"""
    function moments!(state::State)
Calculates presure and velocity
"""
function moments!(state:: State, sysc:: SysConst)
    f0, f1, f2, f3, f4, f5, f6, f7, f8= LBE.viewdists(state.fout)
    sum!(state.rho, state.fout)
    state.velx .= (sysc.vx[1]*f0 .+ sysc.vx[2]*f1 .+ sysc.vx[3]*f2 .+ sysc.vx[4]*f3 .+ sysc.vx[5]*f4 .+ sysc.vx[6]*f5 .+ sysc.vx[7]*f6 .+ sysc.vx[8]*f7 .+ sysc.vx[9]*f8) ./ state.rho
    state.vely .= (sysc.vy[1]*f0 .+ sysc.vy[2]*f1 .+ sysc.vy[3]*f2 .+ sysc.vy[4]*f3 .+ sysc.vy[5]*f4 .+ sysc.vy[6]*f5 .+ sysc.vy[7]*f6 .+ sysc.vy[8]*f7 .+ sysc.vy[9]*f8) ./ state.rho
    return nothing
end

function moments!(state:: StateWithBound, sysc:: SysConstWithBound)
    f0, f1, f2, f3, f4, f5, f6, f7, f8= LBE.viewdists(state.fout)
    sum!(state.rho, state.fout)
    state.velx .= (sysc.vx[1]*f0 .+ sysc.vx[2]*f1 .+ sysc.vx[3]*f2 .+ sysc.vx[4]*f3 .+ sysc.vx[5]*f4 .+ sysc.vx[6]*f5 .+ sysc.vx[7]*f6 .+ sysc.vx[8]*f7 .+ sysc.vx[9]*f8) ./ state.rho
    state.vely .= (sysc.vy[1]*f0 .+ sysc.vy[2]*f1 .+ sysc.vy[3]*f2 .+ sysc.vy[4]*f3 .+ sysc.vy[5]*f4 .+ sysc.vy[6]*f5 .+ sysc.vy[7]*f6 .+ sysc.vy[8]*f7 .+ sysc.vy[9]*f8) ./ state.rho
    return nothing
end