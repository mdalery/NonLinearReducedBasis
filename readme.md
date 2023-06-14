# NonLinearReducedBasis

Before doing anything please ***read until the Constants section.***

To run the simulations, run Julia from the directory *NonLinearReducedBasis*, then  
- `include("scripts/initialize.jl")` to initialize modules and computational constants,  
- `include("scripts/selections_and_errors.jl")` to run the computations and  
- `include("scripts/plots.jl")` to save the plots.

## Constants

Different constants are defined in scripts/initialize.jl that you should probably change before doing anything.
- `const _NTHREADS` number of threads you want to use, please change this value to match your computer characteristics
- `const SAVE_STRING` directory where results are saved
- `const SAVE_STRING_FIGURES` directory where images are saved
- `const SNAPSHOTS` set $`\mathcal{M}_{tr}`$ of solutions snapshots
- `const PREVIOUSLY_COMPUTED_BASIS` `true` to continue already started computations, find the already computed reduced basis in directory `SAVE_STRING`
- `const N_OFFLINE_SELECT` number of elements to add in the reduced basis
- `const ITERATIONS` maximum number of iterations the LBFGS algorithm makes in the online phase
- `const LENGTH_SOBOL` number of starting point for the online phase
- `const R` size of the box, $B$ in the paper
- `const SMOOTH` small parameter to smoothen the energy functional
- `const PREVIOUSLY_COMPUTED_ONLINE_ERRORS` `true` to continue already started computations, find the already computed energy errors in directory `SAVE_STRING`
- `const MOLECULES` online test set of molecules
- `const N_MAX_ONLINE` maximum number of elements to compute online error on
- `_INTER` and `_EXTER` constants are just two other testing sets (taken outside $`\mathcal{M}_{\mathbf{z}}`$ in the paper)

# Included modules

See *NonLinearReducedBasis/examples/* for examples on the two modules. Loading of modules is in *examples/initialize.jl*.

## MixtureBarycenter.jl

Module for computing the mixture barycenters, see *examples/mixtures.jl* for a small example on how to use.
To add a new type of mixture, add a file in *MixtureBarycenter.jl/atoms* and follow what is done in *MixtureBarycenter.jl/atoms/slater.jl*:
- create a new `Atom` subtype, for example `Gaussian <: Atom`
- add its distribution function as `function ( atom::Gaussian )( x::Float64 )`
- add the atomic distance, could be Wasserstein, as `function sq_atomic_distance( atom1::Gaussian, atom2::Gaussian)`
- add the corresponding barycenter as `function barycenter( atoms::Vector{Gaussian}, λ::Vector{Float64} )`, it needs to return a `Gaussian`
- (optional) add constructors for mixtures of your atomic distribution `Gaussian`

## Greedy.jl

Module for computing solutions of a given molecule, selecting a reduced basis of solution via the greedy algorithm, and computing optimal weights for barycenters as approximations.
For examples, see *examples/greedy.jl* for the offline phase and *examples/energy_projection.jl* for the online phase.
