"""
	Molecule
	Molecule( zs::Vector{Float64}, rs::Vector{Float64} )

	Object describing an unidimensional molecule with charges `zs` and positions `rs`.
"""
struct Molecule
	zs::Vector{Float64}
	rs::Vector{Float64}

	function Molecule( zs::Vector{Float64}, rs::Vector{Float64} )
		@assert length(zs) == length(rs)
		@assert !( false in (zs .> zero(Float64)) )

		new(zs, rs)
	end
end


"""
	matC( molecule::Molecule, ζ::Float64 )

	Compute the matrix ``C_ζ``.
"""
@inline function matC( molecule::Molecule, ζ::Float64 )
	return [ z2 * exp(-ζ * abs(r1 - r2)) for r1 in molecule.rs, (z2, r2) in collect(zip(molecule.zs, molecule.rs)) ]
end


"""
	dettC( molecule::Molecule, ζ::Float64 )

	Compute the determinant of the matrix ``C_ζ  - ζI``.
"""
@inline function detC( molecule::Molecule, ζ::Float64 )
	return det(matC(molecule, ζ) - UniformScaling(ζ))
end


"""
	secante( f::Function, x::Float64, y::Float64; ε::Float64 = 1e-16 )

	Find a zero of `f` via a secante algorithm.
"""
@inline function secante( f::Function, x::Float64, y::Float64; ε::Float64 = 1e-16 )
	@assert ε > zero(Float64)

	while abs(x - y) > ε
		x, y = x - f(x) / ((f(x) - f(y)) / (x - y)), x
	end
	return x
end


"""
	solution( molecule::Molecule; ε::Float64 = 1e-16 )

	Compute the ground state of a given `molecule` as a `Mixture{Slater}`.
"""
function solution( molecule::Molecule; ε::Float64 = 1e-16 )
	@assert ε > zero(Float64)

	ζ = secante(ζ -> detC(molecule, ζ), sum(molecule.zs), sum(molecule.zs) - 0.01; ε)
	u = eigvecs(matC(molecule, ζ))[:, end] .* molecule.zs
	if @inbounds u[1] < 0.0
		u = -u
	end
	u = u / sum(u)
	return MixtureSlater(ζ, molecule.rs, u)
end
