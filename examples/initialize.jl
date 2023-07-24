@info "Initalizing packages..."

using Distributed
using JLD2

const _NTHREADS = 1

if nprocs() == 1
	addprocs(_NTHREADS) # Number of threads available
end

@everywhere include("../MixtureBarycenter.jl/MixtureBarycenter.jl")
@everywhere include("../Greedy.jl/Greedy.jl")
using Main.MixtureBarycenter
@everywhere using Main.Greedy

@everywhere GC.gc()

@info "Initalization done."
