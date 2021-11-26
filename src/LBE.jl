module LBE

using Plots

include("initialize.jl")
include("logic.jl")
include("moments.jl")
include("equilibrium.jl")
include("collide.jl")
include("forcing.jl")
include("initialvalues.jl")
include("measure.jl")
include("simulate.jl")

end