using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("initialize.jl")
include("selection_and_errors.jl")
include("plots.jl")
