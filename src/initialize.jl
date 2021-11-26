abstract type LBM_state end

"""
    SyConst{T}

Struct that contains all run time constants

# Arguments

TODO: add Docstring
"""
Base.@kwdef struct SysConst{T}
    Lx:: Int = 256
    Ly:: Int = 256
    #Lz if D = 3
    Tmax:: Int = 1000
    tdump::Int = Tmax/10
    D:: Int = 2
    Q:: Int = 9

    tau:: T = 1.0
    utau:: T = 1.0/tau
    umtau:: T = 1.0- utau
    c_s:: T = 1.0/sqrt(3.0)
    csq:: T = c_s*c_s
    twocpow4:: T = 2* c_s*c_s*c_s*c_s
    twocsq:: T = 2*c_s*c_s
    # Here one could implement value dependenon D and Q at the moment it is D2Q9
    w:: Vector{T} = [4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0]
    vx:: Vector{T} = [0.0,1.0,1.0, 0.0, -1, -1.0, -1, 0.0, 1]
    vy:: Vector{T} = [0.0,0.0,1.0, 1.0, 1, 0.0, -1, -1.0, -1]
end

Base.@kwdef struct SysConstWithBound{T}
    Lx:: Int = 256
    Ly:: Int = 256
    #Lz if D = 3
    Tmax:: Int = 1000
    tdump::Int = Tmax/10
    D:: Int = 2
    Q:: Int = 9

    tau:: T = 1
    utau:: T = 1.0/tau
    umtau:: T = 1.0- utau
    c_s:: T = 1.0/sqrt(3.0)
    csq:: T = c_s*c_s
    twocpow4:: T = 2* c_s*c_s*c_s*c_s
    twocsq:: T = 2*c_s*c_s
    # Here one could implement value dependenon D and Q at the moment it is D2Q9
    w:: Vector{T} = [4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0]
    vx:: Vector{T} = [0.0,1.0,1.0, 0.0, -1, -1.0, -1, 0.0, 1]
    vy:: Vector{T} = [0.0,0.0,1.0, 1.0, 1, 0.0, -1, -1.0, -1]


    obs::Array{T,2} = zeros(Lx,Ly)
    obslist = obslist(obs, T=T)
    interior::Array{T,2} = obslist[1]
    obsleft::Array{T,2} = obslist[2]
    obsright::Array{T,2} =obslist[3]
    obsup::Array{T,2} = obslist[4]
    obsdown::Array{T,2} = obslist[5]
    cornerru::Array{T,2} = obsright .* obsup
    cornerrd::Array{T,2} = obsright .* obsdown
    cornerlu::Array{T,2} = obsleft .* obsup
    cornerld::Array{T,2} = obsleft .* obsdown

end




"""
    State{T,N}
Storing all values of simulation

# Arguments

- `fout :: Array{T,N}`: Output distribution function
- `ftemp :: Array{T,N}`: Temporary distribution function, only used if `sys.τ ≠ 1`
- `feq :: Array{T,N}`: Equilibrium distribution function
- `rho :: Matrix{T}`: Field that stores the scalar density values
- `velx :: Matrix{T}`: Field that stores the x-component of the velocity vector
- `vely :: Matrix{T}`: Field that stores the y-component of the velocity vector
- `vsq :: Matrix{T}`: Field that stores the velocity squared, used in `equilibrium!`
""" 
Base.@kwdef struct State{T, N} <: LBM_state
    fout:: Array{T, N}
    ftemp:: Array{T,N}
    feq:: Array{T,N}
    p:: Array{T,N}

    rho:: Matrix{T}
    velx:: Matrix{T}
    vely:: Matrix{T}
    velsq:: Matrix{T}
end

Base.@kwdef struct StateWithBound{T, N} <: LBM_state
    fout:: Array{T, N}
    ftemp:: Array{T,N}
    feq:: Array{T,N}
    fbound:: Array{T, N}
    fboundinversed:: Array{T,N}
    p:: Array{T,N}

    rho:: Matrix{T}
    velx:: Matrix{T}
    vely:: Matrix{T}
    velsq:: Matrix{T}
end


"""
    function Sys(sysc::SysConst)
space allocation
"""
function MyState(sysc::SysConst; T=Float64)
    state = State{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly),
        velx = zeros(sysc.Lx, sysc.Ly),
        vely = zeros(sysc.Lx, sysc.Ly),
        velsq = zeros(sysc.Lx, sysc.Ly)
    )
    return state
end


function MyState(sysc::SysConstWithBound; T=Float64)
    state = StateWithBound{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fbound= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fboundinversed= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly)-sysc.interior,
        velx = zeros(sysc.Lx, sysc.Ly),
        vely = zeros(sysc.Lx, sysc.Ly),
        velsq = zeros(sysc.Lx, sysc.Ly)
    )
    return state
end

function MyStateX(sysc::SysConst; T=Float64)
    state = State{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly),
        velx = ones(sysc.Lx, sysc.Ly) ,
        vely = zeros(sysc.Lx, sysc.Ly),
        velsq = ones(sysc.Lx, sysc.Ly)
    )
    return state
end

function MyStateY(sysc::SysConst; T=Float64)
    state = State{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly),
        velx = zeros(sysc.Lx, sysc.Ly) ,
        vely = ones(sysc.Lx, sysc.Ly),
        velsq = ones(sysc.Lx, sysc.Ly)
    )
    return state
end



function MyStateX(sysc::SysConstWithBound; T=Float64)
    state = StateWithBound{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fbound= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fboundinversed= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly)-sysc.interior,
        velx = ones(sysc.Lx, sysc.Ly),
        vely = zeros(sysc.Lx, sysc.Ly),
        velsq = ones(sysc.Lx, sysc.Ly)
    )
    return state
end


function MyStateY(sysc::SysConstWithBound; T=Float64)
    state = StateWithBound{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fbound= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fboundinversed= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly)-sysc.interior,
        velx = zeros(sysc.Lx, sysc.Ly),
        vely = ones(sysc.Lx, sysc.Ly),
        velsq = ones(sysc.Lx, sysc.Ly)
    )
    return state
end
"""
    function obslist(obs)
Returns a list containing obsup, obsdown, obsleft, obsright, interior

TODO: Find out how double corners e.g. leftup and leftdown should be handled 
"""
function obslist(obs; T=Float64)
    Lx = size(obs)[1]
    Ly = size(obs)[2]
    interior::Array{T,2} = zeros(Lx,Ly)
    obsleft::Array{T,2} = zeros(Lx, Ly)
    obsright::Array{T,2} = zeros(Lx, Ly)
    obsup::Array{T,2} = zeros(Lx, Ly)
    obsdown::Array{T,2} = zeros(Lx, Ly)
    for i in 1:Lx, j in 1:Ly
        interior[i,j] = (obs[i,j]==1 && obs[i,LBE.previndex(j,Ly)]==1  && obs[i,LBE.nextindex(j,Ly)]==1 && obs[LBE.previndex(i,Lx),j]==1 && obs[LBE.previndex(i,Lx),LBE.previndex(j,Ly)]==1 && obs[LBE.previndex(i,Lx),LBE.nextindex(j,Ly)]==1 && obs[LBE.nextindex(i,Lx),j]==1 && obs[LBE.nextindex(i,Lx),LBE.previndex(j,Ly)]==1 && obs[LBE.nextindex(i,Lx),LBE.nextindex(j,Ly)]==1) ? 1 : 0 
        #-------------------------------------Linesegment-----------------------------------------corner----------------------------------left up corner-------------------------------------------------------------------------------left down corner--------------------------------------------------
        obsleft[i,j] = (obs[i,j]==1 && ((obs[i, LBE.nextindex(j,Ly)]==0) || (obs[i, LBE.nextindex(j,Ly)]==1 && (  (obs[LBE.nextindex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[LBE.nextindex(i,Lx),j]==1)   || (obs[LBE.previndex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[LBE.previndex(i,Lx),j]==1) )))) ? 1 : 0
        #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right up corner-------------------------------------------------------------------------------right down corner--------------------------------------------------
        obsright[i,j] = (obs[i,j]==1 && ((obs[i, LBE.previndex(j,Ly)]==0) || (obs[i, LBE.previndex(j,Ly)]==1 && (  (obs[LBE.nextindex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[LBE.nextindex(i,Lx),j]==1)   || (obs[LBE.previndex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[LBE.previndex(i,Lx),j]==1) )))) ? 1 : 0
        #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right up corner-------------------------------------------------------------------------------left up corner--------------------------------------------------
        obsup[i,j] = (obs[i,j]==1 && ((obs[LBE.nextindex(i,Lx), j]==0) || (obs[LBE.nextindex(i,Lx), j]==1 && (  (obs[LBE.nextindex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[i,LBE.previndex(j,Ly)]==1)   || (obs[LBE.nextindex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[i,LBE.nextindex(j,Ly)]==1)  )))) ? 1 : 0
        #-------------------------------------Linesegment-----------------------------------------corner----------------------------------right down corner-----------------------------------------------------------------------------left down corner--------------------------------------------------
        obsdown[i,j] = (obs[i,j]==1 && ((obs[LBE.previndex(i,Lx), j]==0) || (obs[LBE.previndex(i,Lx), j]==1 && (  (obs[LBE.previndex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[i,LBE.previndex(j,Ly)]==1)   || (obs[LBE.previndex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[i,LBE.nextindex(j,Ly)]==1)  )))) ? 1 : 0
    end
    return [interior, obsleft, obsright, obsup, obsdown]
end

"""
    function initdis!(sys::SysConst, state::State)
for nonzero velocity field you need to initialize the distribution as f_i = f_i^eq. 
    Be carefull, that is not sufficent for incompressible fluids as rho and v are not independent there
"""
function initdis!(state::State, sys::SysConst)
    equilibrium!(state, sys)
    state.ftemp .= state.feq
    state.fout .= state.feq
    return nothing
end


function initdis!(state::StateWithBound, sys::SysConstWithBound)
    equilibrium!(state, sys)
    state.ftemp .= state.feq
    state.fout .= state.feq
    return nothing
end

