"""
    function BGKandStream!(state:: State, sys::SysConst)
Bhatnagar-Gross-Krook equation consiting of collisopn and streaming step
"""
function BGKandStream!(state:: State, sys::SysConst)
    fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8 = viewdists(state.feq)
    ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8 = viewdists(state.ftemp)
    fo0, fo1, fo2, fo3, fo4, fo5, fo6, fo7, fo8 = viewdists(state.fout)

    #colide
    fo0 .= ft0 * sys.umtau .+ fe0*sys.utau
    fo1 .= ft1 * sys.umtau .+ fe1*sys.utau
    fo2 .= ft2 * sys.umtau .+ fe2*sys.utau
    fo3 .= ft3 * sys.umtau .+ fe3*sys.utau
    fo4 .= ft4 * sys.umtau .+ fe4*sys.utau
    fo5 .= ft5 * sys.umtau .+ fe5*sys.utau
    fo6 .= ft6 * sys.umtau .+ fe6*sys.utau
    fo7 .= ft7 * sys.umtau .+ fe7*sys.utau
    fo8 .= ft8 * sys.umtau .+ fe8*sys.utau

    #stream fti becomes the circhift of foi in given direction
    circshift!(ft0, fo0, (0, 0))
    circshift!(ft1, fo1, (0,1))
    circshift!(ft2, fo2, (-1,1))
    circshift!(ft3, fo3, (-1,0))
    circshift!(ft4, fo4, (-1,-1))
    circshift!(ft5, fo5, (0,-1))
    circshift!(ft6, fo6, (1, -1))
    circshift!(ft7, fo7, (1,0))
    circshift!(ft8, fo8, (1,1))

    state.fout .= state.ftemp
    return nothing
end


function BGKandStreamChannel!(state:: StateWithBound, sys::SysConstWithBound)
    fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8 = LBE.viewdists(state.feq)
    ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8 = LBE.viewdists(state.ftemp)
    fo0, fo1, fo2, fo3, fo4, fo5, fo6, fo7, fo8 = LBE.viewdists(state.fout)
    fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7, fb8 = LBE.viewdists(state.fbound)
    fbi0, fbi1, fbi2, fbi3, fbi4, fbi5, fbi6, fbi7, fbi8 = LBE.viewdists(state.fboundinversed)

    

    #colide
    fo0 .= ft0 * sys.umtau .+ fe0*sys.utau
    fo1 .= ft1 * sys.umtau .+ fe1*sys.utau 
    fo2 .= ft2 * sys.umtau .+ fe2*sys.utau 
    fo3 .= ft3 * sys.umtau .+ fe3*sys.utau
    fo4 .= ft4 * sys.umtau .+ fe4*sys.utau
    fo5 .= ft5 * sys.umtau .+ fe5*sys.utau
    fo6 .= ft6 * sys.umtau .+ fe6*sys.utau
    fo7 .= ft7 * sys.umtau .+ fe7*sys.utau
    fo8 .= ft8 * sys.umtau .+ fe8*sys.utau

    #storing uper bounds
    fb2 .= fo2 .* sys.obsup
    fb3 .= fo3 .* sys.obsup
    fb4 .= fo4 .* sys.obsup

    #storing lower bounds
    fb6 .= fo6 .* sys.obsdown
    fb7 .= fo7 .* sys.obsdown
    fb8 .= fo8 .* sys.obsdown


    #storing right bounds
    fb2 .+= fo2 .* (sys.obsright .- sys.cornerru) 
    fb1 .= fo1 .* sys.obsright
    fb8 .+= fo8 .* (sys.obsright .- sys.cornerrd)

    #storing left bounds
    fb4 .+= fo4 .* (sys.obsleft .- sys.cornerlu)
    fb5 .= fo5 .* sys.obsleft
    fb6 .+= fo6 .* (sys.obsleft .- sys.cornerld)

    #erasing the boundary distribution values
    fo1 .-= fb1
    fo2 .-= fb2
    fo3 .-= fb3
    fo4 .-= fb4
    fo5 .-= fb5
    fo6 .-= fb6
    fo7 .-= fb7
    fo8 .-= fb8


    #stream fti becomes the circhift of foi in given direction
    circshift!(ft0, fo0, (0, 0))
    circshift!(ft1, fo1, (0,1))
    circshift!(ft2, fo2, (-1,1))
    circshift!(ft3, fo3, (-1,0))
    circshift!(ft4, fo4, (-1,-1))
    circshift!(ft5, fo5, (0,-1))
    circshift!(ft6, fo6, (1, -1))
    circshift!(ft7, fo7, (1,0))
    circshift!(ft8, fo8, (1,1))

    #interchange stored values
    #up
    fbi6 .= fb2 .* sys.obsup
    fbi7 .= fb3 .* sys.obsup
    fbi8 .= fb4 .* sys.obsup

    #down
    fbi2 .= fb6 .* sys.obsdown
    fbi3 .= fb7 .* sys.obsdown
    fbi4 .= fb8 .* sys.obsdown

    #right
    fbi6 .+= fb2 .* (sys.obsright .- sys.cornerru)
    fbi5 .= fb1 .* sys.obsright
    fbi4 .+= fb8 .* (sys.obsright  .- sys.cornerrd)


    #left
    fbi8 .+= fb4 .* (sys.obsleft .- sys.cornerlu)
    fbi1 .= fb5 .* sys.obsleft
    fbi2 .+= fb6 .* (sys.obsleft .- sys.cornerld)

    #add the stored boundary terms
    ft1 .+= fbi1
    ft2 .+= fbi2
    ft3 .+= fbi3
    ft4 .+= fbi4
    ft5 .+= fbi5
    ft6 .+= fbi6
    ft7 .+= fbi7
    ft8 .+= fbi8






    state.fout .= state.ftemp
    return nothing
end