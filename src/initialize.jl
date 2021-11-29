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
  # The names where given as follows

#
#  For a empty box
#        xxxxxxxx  obsup
#        x      x
# obsleftx      x obsright
#        xxxxxxxx obsdown
#  
# Which makes the names weird for a filled shape
#           obsdown
#          xxxxxxxx  
# obsright xxxxxxxx obsleft
#          xxxxxxxx 
#          xxxxxxxx 
#           obsup
#Corners are named like this
# corneroutlu      corneroutru
#          xxxxxxxx  
#          xxxxxxxx
#          xxxxxxxx 
#          xxxxxxxx 
# corneroutld      corneroutrd
#
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# x cornerlu              cornerru  x 
# x                                 x
# x cornerld               cornerrd x 
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    obs::Array{T,2} = zeros(Lx,Ly)
    obslist = obslist(obs, T=T)
    # interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld
    interior::Array{T,2} = obslist[1]
    border = obslist[2]
    # obsleft::Array{T,2} = obslist[2]
    # obsright::Array{T,2} =obslist[3]
    # obsup::Array{T,2} = obslist[4]
    # obsdown::Array{T,2} = obslist[5]
    # corneroutlu::Array{T,2} = obslist[6]
    # corneroutld::Array{T,2} = obslist[7]
    # corneroutru::Array{T,2} = obslist[8]
    # corneroutrd::Array{T,2} = obslist[9]
    # cornerru::Array{T,2} = obslist[10]
    # cornerrd::Array{T,2} = obslist[11]
    # cornerlu::Array{T,2} = obslist[12]
    # cornerld::Array{T,2} = obslist[13]

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


"""
 function MyStateX(sysc, kwargs...)#
creates velocity field in X direction
TODO: initialize with field and state, so that you can create whatever vel field you want
"""
function MyStateX(sysc::SysConstWithBound; T=Float64)
    state = StateWithBound{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fbound= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        p=zeros(sysc.Lx, sysc.Ly, sysc.Q),

        rho= ones(sysc.Lx, sysc.Ly)-sysc.interior,
        velx = ones(sysc.Lx, sysc.Ly),
        vely = zeros(sysc.Lx, sysc.Ly),
        velsq = ones(sysc.Lx, sysc.Ly)
    )
    return state
end

"""
 function MyStateY(sysc, kwargs...)#
creates velocity field in Y direction
TODO: initialize with field and state, so that you can create whatever vel field you want
"""
function MyStateY(sysc::SysConstWithBound; T=Float64)
    state = StateWithBound{T, 3}(
        fout = zeros(sysc.Lx, sysc.Ly, sysc.Q),
        ftemp= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        feq= zeros(sysc.Lx, sysc.Ly, sysc.Q),
        fbound= zeros(sysc.Lx, sysc.Ly, sysc.Q),
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
Returns a list containing interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld

TODO: Find out how double corners e.g. leftup and leftdown should be handled 
        Add a field for outer corners. 
"""
function obslist(obs; T=Float64)
    Lx = size(obs)[1]
    Ly = size(obs)[2]
    interior::Array{T,2} = zeros(Lx,Ly)
    obsleft::Array{T,2} = zeros(Lx, Ly)
    obsright::Array{T,2} = zeros(Lx, Ly)
    obsup::Array{T,2} = zeros(Lx, Ly)
    obsdown::Array{T,2} = zeros(Lx, Ly)
    corneroutru::Array{T,2} = zeros(Lx, Ly)
    corneroutrd::Array{T,2} = zeros(Lx, Ly)
    corneroutlu::Array{T,2} = zeros(Lx, Ly)
    corneroutld::Array{T,2} = zeros(Lx, Ly)
    cornerru::Array{T,2} = zeros(Lx, Ly)
    cornerrd::Array{T,2} = zeros(Lx, Ly)
    cornerlu::Array{T,2} = zeros(Lx, Ly)
    cornerld::Array{T,2} = zeros(Lx, Ly)

      # The names where given as follows

#
#  For a empty box
#        xxxxxxxx  obsup
#        x      x
# obsleftx      x obsright
#        xxxxxxxx obsdown
#  
# Which makes the names weird for a filled shape
#           obsdown
#          xxxxxxxx  
# obsright xxxxxxxx obsleft
#          xxxxxxxx 
#          xxxxxxxx 
#           obsup
#Corners are named like this
# corneroutlu      corneroutru
#          xxxxxxxx  
#          xxxxxxxx
#          xxxxxxxx 
#          xxxxxxxx 
# corneroutld      corneroutrd
#
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# x cornerlu              cornerru  x 
# x                                 X
# x cornerld               cornerrd x 
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
        #outer corners
        corneroutlu[i,j] = ( obs[i,j]==1 && obs[LBE.previndex(i,Lx),j]==0 &&  obs[LBE.previndex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[i,LBE.previndex(j,Ly)]==0) ? 1 : 0
        corneroutld[i,j] = ( obs[i,j]==1 && obs[LBE.nextindex(i,Lx),j]==0 &&  obs[LBE.nextindex(i,Lx),LBE.previndex(j,Ly)]==0 && obs[i,LBE.previndex(j,Ly)]==0) ? 1 : 0
        corneroutru[i,j] = ( obs[i,j]==1 && obs[LBE.previndex(i,Lx),j]==0 &&  obs[LBE.previndex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[i,LBE.nextindex(j,Ly)]==0) ? 1 : 0
        corneroutrd[i,j] = ( obs[i,j]==1 && obs[LBE.nextindex(i,Lx),j]==0 &&  obs[LBE.nextindex(i,Lx),LBE.nextindex(j,Ly)]==0 && obs[i,LBE.nextindex(j,Ly)]==0) ? 1 : 0
    end
    #inner corners
    cornerlu = obsleft .* obsup .- corneroutrd
    cornerld = obsleft .* obsdown .- corneroutru
    cornerru = obsright .* obsup .- corneroutld
    cornerrd = obsright .* obsdown .- corneroutlu

    
    #preparing return
    border1 = zeros(Lx, Ly)
    border2 = zeros(Lx, Ly)
    border3 = zeros(Lx, Ly)
    border4 = zeros(Lx, Ly)
    border5 = zeros(Lx, Ly)
    border6 = zeros(Lx, Ly)
    border7 = zeros(Lx, Ly)
    border8 = zeros(Lx, Ly)     

    border1 .= obsright .- corneroutld .- corneroutlu
    border2 .= obsright .+ obsup .- cornerru .- corneroutld .- corneroutrd .- corneroutlu
    border3 .= obsup .- corneroutld .- corneroutrd
    border4 .= obsleft .+ obsup .- cornerlu .- corneroutld .- corneroutrd .- corneroutru
    border5 .= obsleft .- corneroutrd .- corneroutru
    border6 .= obsleft .+ obsdown .- cornerld .- corneroutlu .- corneroutrd .- corneroutru
    border7 .= obsdown .- corneroutru .- corneroutlu
    border8 .= obsright .+ obsdown .- cornerrd .- corneroutld .- corneroutlu .- corneroutru

     #retrun  interior, obsleft, obsright, obsup, obsdown, corneroutlu, corneroutld, corneroutru, corneroutrd, cornerru, cornerrd, cornerlu, cornerld, [border1, border2, border3, border4, border5, border6, border7, border8] 
    return interior, [border1, border2, border3, border4, border5, border6, border7, border8] 
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

