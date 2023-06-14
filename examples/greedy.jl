## GREEDY SELECTION EXAMPLE

SAVE_STRING_EXAMPLE = "examples/greed_example"

# Compute a solution
Zs = [ 0.8, 1.1 ]
molecule = Molecule(Zs, [ -1.0, 1.0 ])
sol = solution(molecule)

# Snapshots initialization
L = 25
snapshots = [ solution(Molecule(Zs, [ -1.5 + 3*rand(), -1.5 + 3*rand() ])) for n in 1:L ] # 25 random snapshots

# Greed 5 elements
reduced = Mixture{Slater}[]
errs_max = Float64[]
errs_mean = Float64[]
greed!(reduced, errs_max, errs_mean, 5, snapshots, SAVE_STRING_EXAMPLE)
xs = -3.0:0.01:3.0

# Projection of the solution on the selected basis.
w = weights_min(reduced)
A = matrixA(reduced, w)
invA = inv(A)
_, λ, _ = projection(sol, reduced, w, invA)
projected5 = barycenter(reduced, λ, w) # Mixture of Slater distribution containing the projection

# Go up to 8 elements
greed!(reduced, errs_max, errs_mean, 3, snapshots, SAVE_STRING_EXAMPLE)

w = weights_min(reduced)
A = matrixA(reduced, w)
invA = inv(A)
cost, λ = projection(sol, reduced, w, invA)
projected8 = barycenter(reduced, λ, w) # Mixture of Slater distribution containing the projection
