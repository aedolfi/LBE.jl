function viewdists(f)
    f0 = view(f, :, :, 1)
    f1 = view(f, :, :, 2)
    f2 = view(f, :, :, 3)
    f3 = view(f, :, :, 4)
    f4 = view(f, :, :, 5)
    f5 = view(f, :, :, 6)
    f6 = view(f, :, :, 7)
    f7 = view(f, :, :, 8)
    f8 = view(f, :, :, 9)
    
    return f0, f1, f2, f3, f4, f5, f6, f7, f8
end



"""
    function previndex(i, L)
gives the index before i on a 1:L torus
"""
function previndex(i,L)
    if i==1
        return L
    else
        return i-1
    end
end


"""
    function nextindex(i, L)
gives the index before i on a 1:L torus
"""
function nextindex(i,L)
    if i==L
        return 1
    else
        return i+1
    end
end