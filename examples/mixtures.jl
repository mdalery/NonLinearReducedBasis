## MIXTURES BARYCENTERS EXAMPLES

# Three mixtures of Slater distributions
m1 = MixtureSlater(1.0, [ -1.5, 0.5 ], [ 0.4, 0.6 ]) # First mixture
m2 = MixtureSlater([ 0.9, 1.1 ], [ -0.5, 1.0 ], [ 0.3, 0.7 ]) # Second mixture
m3 = MixtureSlater(0.75 + 0.5*rand(), [ -1.5 + 3*rand(), -1.5 + 3*rand() ], [ 0.5, 0.5 ]) # Third mixture (with a bit of randomness)

mixtures = [ m1, m2, m3 ]
w = weights_min(mixtures) # Computation of the weights `w^*` for barycenters

λ = [ 0.3, 0.5, 0.4 ]
b = barycenter(mixtures, λ, w) # Computation of a barycenter between the three above mixtures

@show b(0.25) # Compute the distribution function of b at 0.25
