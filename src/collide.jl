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

#TODO: Fix propagation for outer corners
function BGKandStreamChannel!(state:: StateWithBound, sys::SysConstWithBound)
    fe0, fe1, fe2, fe3, fe4, fe5, fe6, fe7, fe8 = LBE.viewdists(state.feq)
    ft0, ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8 = LBE.viewdists(state.ftemp)
    fo0, fo1, fo2, fo3, fo4, fo5, fo6, fo7, fo8 = LBE.viewdists(state.fout)
    fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7, fb8 = LBE.viewdists(state.fbound)

    

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


        # The names where given as follows

#
#  For a empty box
#        xxxxxxxx  obsup
#        x      x
# obsleftx      x obsright
#        x      x
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


#TODO: In order to optimize the algorithm all the following differences should be saved in SysConstWithBound

    #Storing bounds that shoulb be handeld differnetly due to corners
    # fb1 .= fo1 .* (sys.obsright .- sys.corneroutld .- sys.corneroutlu)
    # fb2 .= fo2 .* (sys.obsright .+ sys.obsup .- sys.cornerru .-sys.corneroutld .- sys.corneroutrd .- sys.corneroutlu) 
    # fb3 .= fo3 .* (sys.obsup .- sys.corneroutld .- sys.corneroutrd)
    # fb4 .= fo4 .* (sys.obsleft .+ sys.obsup .- sys.cornerlu .- sys.corneroutld .- sys.corneroutrd .- sys.corneroutru)
    # fb5 .= fo5 .* (sys.obsleft .- sys.corneroutrd .- sys.corneroutru)
    # fb6 .= fo6 .* (sys.obsleft .+ sys.obsdown .- sys.cornerld .- sys.corneroutlu .- sys.corneroutrd .- sys.corneroutru)
    # fb7 .= fo7 .* (sys.obsdown .- sys.corneroutru .- sys.corneroutlu)
    # fb8 .= fo8 .* (sys.obsright .+ sys.obsdown .- sys.cornerrd .- sys.corneroutld .- sys.corneroutlu .- sys.corneroutru)


    fb1 .= fo1 .* sys.border[1]
    fb2 .= fo2 .* sys.border[2]
    fb3 .= fo3 .* sys.border[3]
    fb4 .= fo4 .* sys.border[4]
    fb5 .= fo5 .* sys.border[5]
    fb6 .= fo6 .* sys.border[6]
    fb7 .= fo7 .* sys.border[7]
    fb8 .= fo8 .* sys.border[8]

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



    #add the stored boundary terms
    ft1 .+= fb5
    ft2 .+= fb6
    ft3 .+= fb7
    ft4 .+= fb8
    ft5 .+= fb1
    ft6 .+= fb2
    ft7 .+= fb3
    ft8 .+= fb4






    state.fout .= state.ftemp
    return nothing
end