"""
    function equilibrium!(state::State, sysc:: SysConst)
computes f^equilibrium_i = w_i rho ( 1 = 3(vel* v_i) + 9/2(vel* v_i) -3/2 vel^2) to solve incompressible euler Equation
"""
function equilibrium!(state::State, sysc::SysConst)
    f0,f1,f2,f3,f4,f5,f6,f7,f8 = LBE.viewdists(state.feq)
    state.velsq .= state.velx .* state.velx .+ state.vely .* state.vely
    p0,p1,p2,p3,p4,p5,p6,p7,p8  = LBE.viewdists(state.p)
    p0 .= sysc.vx[1] .* state.velx .+ sysc.vy[1] .* state.vely
    p1 .= sysc.vx[2] .* state.velx .+ sysc.vy[2] .* state.vely
    p2 .= sysc.vx[3] .* state.velx .+ sysc.vy[3] .* state.vely
    p3 .= sysc.vx[4] .* state.velx .+ sysc.vy[4] .* state.vely
    p4 .= sysc.vx[5] .* state.velx .+ sysc.vy[5] .* state.vely
    p5 .= sysc.vx[6] .* state.velx .+ sysc.vy[6] .* state.vely
    p6 .= sysc.vx[7] .* state.velx .+ sysc.vy[7] .* state.vely
    p7 .= sysc.vx[8] .* state.velx .+ sysc.vy[8] .* state.vely
    p8 .= sysc.vx[9] .* state.velx .+ sysc.vy[9] .* state.vely
    f0 .= sysc.w[1] .* state.rho .*(1.0 .- 1.5 .* state.velsq)
    f1 .= sysc.w[2] .* state.rho .*(1.0 .+ 3.0 .* p1 .+ 4.5 .* p1 .* p1 .- 1.5 .* state.velsq)
    f2 .= sysc.w[3] .* state.rho .*(1.0 .+ 3.0 .* p2 .+ 4.5 .* p2 .* p2 .- 1.5 .* state.velsq)
    f3 .= sysc.w[4] .* state.rho .*(1.0 .+ 3.0 .* p3 .+ 4.5 .* p3 .* p3 .- 1.5 .* state.velsq)
    f4 .= sysc.w[5] .* state.rho .*(1.0 .+ 3.0 .* p4 .+ 4.5 .* p4 .* p4 .- 1.5 .* state.velsq)
    f5 .= sysc.w[6] .* state.rho .*(1.0 .+ 3.0 .* p5 .+ 4.5 .* p5 .* p5 .- 1.5 .* state.velsq)
    f6 .= sysc.w[7] .* state.rho .*(1.0 .+ 3.0 .* p6 .+ 4.5 .* p6 .* p6 .- 1.5 .* state.velsq)
    f7 .= sysc.w[8] .* state.rho .*(1.0 .+ 3.0 .* p7 .+ 4.5 .* p7 .* p7 .- 1.5 .* state.velsq)
    f8 .= sysc.w[9] .* state.rho .*(1.0 .+ 3.0 .* p8 .+ 4.5 .* p8 .* p8 .- 1.5 .* state.velsq)
    return nothing
end

function equilibrium!(state::StateWithBound, sysc::SysConstWithBound)
    f0,f1,f2,f3,f4,f5,f6,f7,f8 = LBE.viewdists(state.feq)
    state.velsq .= state.velx .* state.velx .+ state.vely .* state.vely
    p0,p1,p2,p3,p4,p5,p6,p7,p8  = LBE.viewdists(state.p)
    p0 .= sysc.vx[1] .* state.velx .+ sysc.vy[1] .* state.vely
    p1 .= sysc.vx[2] .* state.velx .+ sysc.vy[2] .* state.vely
    p2 .= sysc.vx[3] .* state.velx .+ sysc.vy[3] .* state.vely
    p3 .= sysc.vx[4] .* state.velx .+ sysc.vy[4] .* state.vely
    p4 .= sysc.vx[5] .* state.velx .+ sysc.vy[5] .* state.vely
    p5 .= sysc.vx[6] .* state.velx .+ sysc.vy[6] .* state.vely
    p6 .= sysc.vx[7] .* state.velx .+ sysc.vy[7] .* state.vely
    p7 .= sysc.vx[8] .* state.velx .+ sysc.vy[8] .* state.vely
    p8 .= sysc.vx[9] .* state.velx .+ sysc.vy[9] .* state.vely
    f0 .= sysc.w[1] .* state.rho .*(1.0 .- 1.5 .* state.velsq)
    f1 .= sysc.w[2] .* state.rho .*(1.0 .+ 3.0 .* p1 .+ 4.5 .* p1 .* p1 .- 1.5 .* state.velsq)
    f2 .= sysc.w[3] .* state.rho .*(1.0 .+ 3.0 .* p2 .+ 4.5 .* p2 .* p2 .- 1.5 .* state.velsq)
    f3 .= sysc.w[4] .* state.rho .*(1.0 .+ 3.0 .* p3 .+ 4.5 .* p3 .* p3 .- 1.5 .* state.velsq)
    f4 .= sysc.w[5] .* state.rho .*(1.0 .+ 3.0 .* p4 .+ 4.5 .* p4 .* p4 .- 1.5 .* state.velsq)
    f5 .= sysc.w[6] .* state.rho .*(1.0 .+ 3.0 .* p5 .+ 4.5 .* p5 .* p5 .- 1.5 .* state.velsq)
    f6 .= sysc.w[7] .* state.rho .*(1.0 .+ 3.0 .* p6 .+ 4.5 .* p6 .* p6 .- 1.5 .* state.velsq)
    f7 .= sysc.w[8] .* state.rho .*(1.0 .+ 3.0 .* p7 .+ 4.5 .* p7 .* p7 .- 1.5 .* state.velsq)
    f8 .= sysc.w[9] .* state.rho .*(1.0 .+ 3.0 .* p8 .+ 4.5 .* p8 .* p8 .- 1.5 .* state.velsq)
    return nothing
end