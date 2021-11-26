# whatever you want to meassure, put your function here


"""
    function vectorfield2d(field, points, arrowlength=0.1)
Plots a vectorfield. 
from: https://discourse.julialang.org/t/plotting-vector-fields/32704
"""
function vectorfield2d(field, points, arrowlength=0.1)
    # More input pattern parsing is solved by the Plots package, but I don't know how.
    errormessage = "Incorrect formatting of points. Please format them as [x1 y1; x2, y2;...]"
    
    if typeof(points) <: Array{<:Number, 2} && size(points)[1] === 2
        vectors = similar(points)
        for i in 1:size(points)[2]
            vectors[:, i] .= collect(field(points[:, i]...))
        end
    else
        error(errormessage)
    end
    vectors .*= arrowlength
    quiver(points[1, :],points[2, :],quiver=(vectors[1, :], vectors[2, :]))
    display(plot!())
end