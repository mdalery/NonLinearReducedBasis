## ENERGY PROJECTION EXAMPLE

# A few solutions
mixtures = [ solution(Molecule([ 0.8, 1.1 ], [ -r, r ])) for r in [ 0.5, 2.5, 1.26 ] ]
w = weights_min(mixtures) # the weights `w^*`

# A molecule whose solution approximation is wanted
molecule = Molecule([ 0.8, 1.1 ], [ -1.78, 1.78 ])

# The barycenter obtain by energy minimization
R = 2.0
len_sobol = 2_000
iterations = 10_000
_, λ = energy_minimization(molecule, mixtures, w, R, len_sobol, iterations)
proj = barycenter(mixtures, λ, w) # The barycenter projection
