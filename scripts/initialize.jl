@info "Initalizing packages and constants..."

using Distributed
using JLD2
using Plots
pyplot()

# Distributed constants
const _NTHREADS = 15 # Number of threads available

if nprocs() == 1
	addprocs(_NTHREADS)
end

@everywhere include("../MixtureBarycenter.jl/MixtureBarycenter.jl")
@everywhere include("../Greedy.jl/Greedy.jl")
using Main.MixtureBarycenter
@everywhere using Main.Greedy

@everywhere GC.gc()

# Saving constants
const SAVE_STRING = "saves/"
const SAVE_STRING_FIGURES = "images/"

# Manifolds constant. If PREVIOUSLY_COMPUTED is true, make sure the below set is the same used before.
const SNAPSHOTS = [ solution(Molecule([0.8, 1.1], [-r, r])) for r in 0.5:0.01:3.0 ]

# Offline constants
const PREVIOUSLY_COMPUTED_BASIS = true
const N_OFFLINE_SELECT = 0 # Put to zero to not compute more

# Online constants
const PREVIOUSLY_COMPUTED_ONLINE_ERRORS = true
const ITERATIONS = 10_000
const LENGTH_SOBOL = 2_000
const R = 2.0
const SMOOTH = 1e-8
const MOLECULES = [ Molecule([ 0.8, 1.1 ], [ -r, r ]) for r in 0.5:0.05:3.0 ]
const N_MAX_ONLINE = 0 # Put to zero to not compute more

# Online internal extrapolation constants
const PREVIOUSLY_COMPUTED_INTER_ERRORS = true
const MOLECULES_INTER = [ [ Molecule([1.9], [0.0]) ]; [ Molecule([0.8, 1.1], [-r, r]) for r in 0.03:0.03:0.48 ] ]
const N_MAX_INTER = 0 # Put to zero to not compute more

# Online external extrapolation constants
const PREVIOUSLY_COMPUTED_EXTER_ERRORS = true
const MOLECULES_EXTER = [ Molecule([0.8, 1.1], [-r, r]) for r in 3.05:0.05:4.0 ]
const N_MAX_EXTER = 0 # Put to zero to not compute more


basis = Mixture{Slater}[]
worsts_proj, means_proj = Float64[], Float64[]
if PREVIOUSLY_COMPUTED_BASIS
	@info "Loading previously computed basis and errors..."
	basis = load_object(SAVE_STRING * "basis.jld2")
	worsts_proj, means_proj = load_object(SAVE_STRING * "projection_errors.jld2")
end

worsts_energy, means_energy = Float64[], Float64[]
if PREVIOUSLY_COMPUTED_ONLINE_ERRORS
	@info "Loading previously computed online errors..."
	worsts_energy, means_energy = load_object(SAVE_STRING * "energy_errors.jld2")
end

worsts_energy_inter, means_energy_inter = Float64[], Float64[]
if PREVIOUSLY_COMPUTED_INTER_ERRORS
	@info "Loading previously computed internal errors..."
	worsts_energy_inter, means_energy_inter = load_object(SAVE_STRING * "inter_energy_errors.jld2")
end

worsts_energy_exter, means_energy_exter = Float64[], Float64[]
if PREVIOUSLY_COMPUTED_EXTER_ERRORS
	@info "Loading previously computed external errors..."
	worsts_energy_exter, means_energy_exter = load_object(SAVE_STRING * "exter_energy_errors.jld2")
end


@info "Initalization done."
