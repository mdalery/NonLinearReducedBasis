"""
	energy( molecule::Molecule, mixture::Mixture{Slater} )

	Compute the energy of the molecular system `molecule` at function point `mixture`.
"""
function energy( molecule::Molecule, mixture::Mixture{Slater} )
	K = length(mixture.weights)
	M = length(molecule.rs)

	ζ = @inbounds mixture.atoms[1].ζ
	E = sum( @inbounds mixture.weights[k] * mixture.weights[l] * (1 - ζ*abs(mixture.atoms[k].r - mixture.atoms[l].r)) * exp(-ζ*abs(mixture.atoms[k].r - mixture.atoms[l].r))
			for k in 1:K, l in 1:K )
	E = ζ^3/8 * E - sum( @inbounds molecule.zs[m] * mixture(molecule.rs[m])^2 for m in 1:M )
	norm2 = ζ/4 * sum( @inbounds mixture.weights[k] * mixture.weights[l] * (1 + ζ*abs(mixture.atoms[k].r - mixture.atoms[l].r)) * exp(-ζ*abs(mixture.atoms[k].r - mixture.atoms[l].r))
					  for k in 1:K, l in 1:K )
	E = E / norm2
	return E
end


"""
	energy( molecule::Molecule, mixture::Mixture{Slater} )

	Compute the energy of the molecular system `molecule` at function point the barycenter between the mixtures `mixtures` with weights `λ`.
	Weights ``w^*`` needs to be previously computed.
"""
function energy( λ::Vector{Float64}, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64} ) where N
	begin
		@assert N isa Int64
		@assert length(λ) == length(mixtures) == N
	end
	M = length(molecule.rs)

	ζ = 1 / sum( @inbounds λ[j] / mixtures[j].atoms[1].ζ for j in 1:N )
	@assert ζ > zero(Float64)

	r( kbar::CartesianIndex{N} ) = sum( @inbounds λ[j] * mixtures[j].atoms[kbar[j]].r for j in 1:N )

	G = sum( wkbar * wlbar * (1 - ζ*abs(r(kbar) - r(lbar))) * exp(-ζ*abs(r(kbar) - r(lbar))) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	S = sum( @inbounds molecule.zs[m] * sum( wkbar * exp(-ζ*abs(r(kbar) - molecule.rs[m])) for (kbar, wkbar) in w )^2 for m in 1:M )
	Num = ζ^2/2 * G - ζ * S
	D = sum( wkbar * wlbar * (1 + ζ*abs(r(kbar) - r(lbar))) * exp(-ζ*abs(r(kbar) - r(lbar))) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	E = Num / D
end


"""
	di_energy( λ::Vector{Float64}, i::Int64, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64} ) where N

	Compute the partial derivative over the variable ``λ_i`` of the energy of the molecular system `molecule`
	at function point the barycenter between the mixtures `mixtures` with weights `λ`. Weights ``w^*`` needs to be previously computed.
"""
function di_energy( λ::Vector{Float64}, i::Int64, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64} ) where N
	begin
		@assert N isa Int64
		@assert length(λ) == length(mixtures) == N
		@assert zero(Int64) <= i <= N
	end
	M = length(molecule.rs)

	ζ = 1 / sum( @inbounds λ[j] / mixtures[j].atoms[1].ζ for j in 1:N )
	@assert ζ > zero(Float64)
	di_ζ = @inbounds - ζ^2 / mixtures[i].atoms[1].ζ

	r( kbar::CartesianIndex{N} ) = sum( @inbounds λ[j] * mixtures[j].atoms[kbar[j]].r for j in 1:N )
	di_r( kbar::CartesianIndex{N} ) = @inbounds mixtures[i].atoms[kbar[i]].r
	
	ζdiff( kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = ζ * abs(r(kbar) - r(lbar))
	di_ζdiff( kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = di_ζ * abs(r(kbar) - r(lbar)) + ζ * (di_r(kbar) - di_r(lbar)) * sign(r(kbar) - r(lbar))
	ζdiffmol( kbar::CartesianIndex{N}, m::Int64 ) = @inbounds ζ * abs(r(kbar) - molecule.rs[m])
	di_ζdiffmol( kbar::CartesianIndex{N}, m::Int64 ) = @inbounds di_ζ * abs(r(kbar) - molecule.rs[m]) + ζ * di_r(kbar) * sign(r(kbar) - molecule.rs[m])

	D = sum( wkbar * wlbar * (1 + ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	di_D = - sum( wkbar * wlbar * di_ζdiff(kbar, lbar) * ζdiff(kbar, lbar) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

	G = sum( wkbar * wlbar * (1 - ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	di_G = - sum( wkbar * wlbar * di_ζdiff(kbar, lbar) * (2 - ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

	S = sum( @inbounds molecule.zs[m] * sum( wkbar * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w )^2 for m in 1:M )
	di_S = -2 * sum( @inbounds molecule.zs[m] * sum( wkbar * di_ζdiffmol(kbar, m) * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w ) * sum( wkbar * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w ) for m in 1:M )

	Num = ζ^2/2 * G - ζ * S
	di_Num = di_ζ*ζ * G + ζ^2/2 * di_G - di_ζ * S - ζ * di_S

	return ( di_Num * D - Num * di_D ) / D^2
end


"""
	senergy( molecule::Molecule, mixture::Mixture{Slater}; smooth::Float64 = 1e-8 )

	Compute a smooth version of the energy of the molecular system `molecule` at function point the barycenter between the mixtures `mixtures` with weights `λ`.
	Weights ``w^*`` needs to be previously computed.
"""
function senergy( λ::Vector{Float64}, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64}; smooth::Float64 = 1e-8 ) where N
	begin
		@assert N isa Int64
		@assert length(λ) == length(mixtures) == N
	end
	M = length(molecule.rs)

	ζ = 1 / sum( @inbounds λ[j] / mixtures[j].atoms[1].ζ for j in 1:N )
	@assert ζ > zero(Float64)

	sabs( x::Float64 ) = smooth_abs(x, smooth = smooth)

	r( kbar::CartesianIndex{N} ) = sum( @inbounds λ[j] * mixtures[j].atoms[kbar[j]].r for j in 1:N )

	G = sum( wkbar * wlbar * (1 - ζ*sabs(r(kbar) - r(lbar))) * exp(-ζ*sabs(r(kbar) - r(lbar))) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	S = sum( @inbounds molecule.zs[m] * sum( wkbar * exp(-ζ*sabs(r(kbar) - molecule.rs[m])) for (kbar, wkbar) in w )^2 for m in 1:M )
	Num = ζ^2/2 * G - ζ * S
	D = sum( wkbar * wlbar * (1 + ζ*sabs(r(kbar) - r(lbar))) * exp(-ζ*sabs(r(kbar) - r(lbar))) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	E = Num / D
end


"""
	di_senergy( λ::Vector{Float64}, i::Int64, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64}; smooth::Float64 = 1e-8 ) where N

	Compute the partial derivative over the variable ``λ_i`` of the smoothed version of the energy of the molecular system `molecule`
	at function point the barycenter between the mixtures `mixtures` with weights `λ`. Weights ``w^*`` needs to be previously computed.
"""
function di_senergy( λ::Vector{Float64}, i::Int64, molecule::Molecule, mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64}; smooth::Float64 = 1e-8 ) where N
	begin
		@assert N isa Int64
		@assert length(λ) == length(mixtures) == N
		@assert zero(Int64) <= i <= N
	end
	M = length(molecule.rs)

	ζ = 1 / sum( @inbounds λ[j] / mixtures[j].atoms[1].ζ for j in 1:N )
	@assert ζ > zero(Float64)
	di_ζ = @inbounds - ζ^2 / mixtures[i].atoms[1].ζ

	r( kbar::CartesianIndex{N} ) = sum( @inbounds λ[j] * mixtures[j].atoms[kbar[j]].r for j in 1:N )
	di_r( kbar::CartesianIndex{N} ) = @inbounds mixtures[i].atoms[kbar[i]].r

	sabs( x::Float64 ) = smooth_abs(x, smooth = smooth)
	ssign( x::Float64 ) = smooth_sign(x, smooth = smooth)
	
	ζdiff( kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = ζ * sabs(r(kbar) - r(lbar))
	di_ζdiff( kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = di_ζ * sabs(r(kbar) - r(lbar)) + ζ * (di_r(kbar) - di_r(lbar)) * ssign(r(kbar) - r(lbar))
	ζdiffmol( kbar::CartesianIndex{N}, m::Int64 ) = @inbounds ζ * sabs(r(kbar) - molecule.rs[m])
	di_ζdiffmol( kbar::CartesianIndex{N}, m::Int64 ) = @inbounds di_ζ * sabs(r(kbar) - molecule.rs[m]) + ζ * di_r(kbar) * ssign(r(kbar) - molecule.rs[m])

	D = sum( wkbar * wlbar * (1 + ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	di_D = - sum( wkbar * wlbar * di_ζdiff(kbar, lbar) * ζdiff(kbar, lbar) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

	G = sum( wkbar * wlbar * (1 - ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
	di_G = - sum( wkbar * wlbar * di_ζdiff(kbar, lbar) * (2 - ζdiff(kbar, lbar)) * exp(-ζdiff(kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

	S = sum( @inbounds molecule.zs[m] * sum( wkbar * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w )^2 for m in 1:M )
	di_S = -2 * sum( @inbounds molecule.zs[m] * sum( wkbar * di_ζdiffmol(kbar, m) * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w ) * sum( wkbar * exp(-ζdiffmol(kbar, m)) for (kbar, wkbar) in w ) for m in 1:M )

	Num = ζ^2/2 * G - ζ * S
	di_Num = di_ζ*ζ * G + ζ^2/2 * di_G - di_ζ * S - ζ * di_S

	return ( di_Num * D - Num * di_D ) / D^2
end


"""
	_step_energy_minimization( λ::Vector{Float64}, basis::Vector{Mixture{Slater}}, f::Function, g!::Function, iterations::Int64; clean::Bool = false )

	Perform an LBFGS algorithm on function `f` of derivative `g!` at starting poitn `λ`.
"""
function _step_energy_minimization( λ::Vector{Float64}, basis::Vector{Mixture{Slater}}, f::Function, g!::Function, iterations::Int64; clean::Bool = false )
	res = optimize(f, g!, λ, LBFGS(), Optim.Options(iterations = iterations))
	cost = Optim.minimum(res)
	minimizer = Optim.minimizer(res)

	if clean
		GC.gc();
	end
	return (cost, minimizer)
end


"""
	energy_minimization( molecule::Molecule, basis::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64}, R::Float64, len_sobol::Int64, iterations::Int64; smooth::Float64 = 1e-10 ) where N

	Find the global minimum of the energy function `senergy`. Weights ``w`` needs to be computed.
"""
function energy_minimization( molecule::Molecule, basis::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64}, R::Float64, len_sobol::Int64, iterations::Int64; smooth::Float64 = 1e-8 ) where N
	begin
		@assert N isa Int64
		@assert length(basis) == N
		@assert R > zero(Float64)
		@assert len_sobol >= 10
		@assert iterations >= 10
	end

	M = length(molecule.rs)
	sabs( x::Float64 ) = smooth_abs(x, smooth = smooth)
	ssign( x::Float64 ) = smooth_sign(x, smooth = smooth)
	r( λ::Vector{Float64}, kbar::CartesianIndex{N} ) = sum( @inbounds λ[j] * basis[j].atoms[kbar[j]].r for j in 1:N )
	di_r( i::Int64, kbar::CartesianIndex{N} ) = @inbounds basis[i].atoms[kbar[i]].r
	ζdiff( λ::Vector{Float64}, ζ::Float64, kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = ζ * sabs(r(λ, kbar) - r(λ, lbar))
	di_ζdiff( λ::Vector{Float64}, ζ::Float64, di_ζ::Float64, i::Int64, kbar::CartesianIndex{N}, lbar::CartesianIndex{N} ) = di_ζ * sabs(r(λ, kbar) - r(λ, lbar)) + ζ * (di_r(i, kbar) - di_r(i, lbar)) * ssign(r(λ, kbar) - r(λ, lbar))
	ζdiffmol( λ::Vector{Float64}, ζ::Float64, kbar::CartesianIndex{N}, m::Int64 ) = ζ * sabs(r(λ, kbar) - molecule.rs[m])
	di_ζdiffmol( λ::Vector{Float64}, ζ::Float64, di_ζ::Float64, i::Int64, kbar::CartesianIndex{N}, m::Int64 ) = di_ζ * sabs(r(λ, kbar) - molecule.rs[m]) + ζ * di_r(i, kbar) * ssign(r(λ, kbar) - molecule.rs[m])

	function f( λ::Vector{Float64} )
		ζ = 1 / sum( @inbounds λ[n] / basis[n].atoms[1].ζ for n in 1:N )
		if ζ <= zero(Float64)
			return 1e7 + sum( @inbounds λ[n] / basis[n].atoms[1].ζ for n in 1:N )^10
		else
			E = 0.5 * ζ^2 * sum( wkbar * wlbar * (1 - ζdiff(λ, ζ, kbar, lbar)) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
			E = E - ζ * sum( molecule.zs[m] * sum( wkbar * exp(-ζdiffmol(λ, ζ, kbar, m)) for (kbar, wkbar) in w )^2 for m in 1:M )
			D = sum( wkbar * wlbar * (1 + ζdiff(λ, ζ, kbar, lbar)) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
			E = E / D
			return E
		end
	end

	function g!( grad::Vector{Float64}, λ::Vector{Float64} )
		ζ = 1 / sum( @inbounds λ[n] / basis[n].atoms[1].ζ for n in 1:N )
		if ζ <= zero(Float64)
			for i in 1:N
				grad[i] = 10.0 / basis[i].atoms[1].ζ * sum( @inbounds λ[n] / basis[n].atoms[1].ζ for n in 1:N )^9
			end
		else
			for i in 1:N
				di_ζ = - ζ^2 / basis[i].atoms[1].ζ

				D = sum( wkbar * wlbar * (1.0 + ζdiff(λ, ζ, kbar, lbar)) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
				di_D = - sum( wkbar * wlbar * di_ζdiff(λ, ζ, di_ζ, i, kbar, lbar) * ζdiff(λ, ζ, kbar, lbar) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

				G = sum( wkbar * wlbar * (1.0 - ζdiff(λ, ζ, kbar, lbar)) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )
				di_G = - sum( wkbar * wlbar * di_ζdiff(λ, ζ, di_ζ, i, kbar, lbar) * (2.0 - ζdiff(λ, ζ, kbar, lbar)) * exp(-ζdiff(λ, ζ, kbar, lbar)) for (kbar, wkbar) in w, (lbar, wlbar) in w )

				S = sum( molecule.zs[m] * sum( wkbar * exp(-ζdiffmol(λ, ζ, kbar, m)) for (kbar, wkbar) in w )^2 for m in 1:M )
				di_S = -2.0 * sum( molecule.zs[m] * sum( wkbar * di_ζdiffmol(λ, ζ, di_ζ, i, kbar, m) * exp(-ζdiffmol(λ, ζ, kbar, m)) for (kbar, wkbar) in w ) * sum( wkbar * exp(-ζdiffmol(λ, ζ, kbar, m)) for (kbar, wkbar) in w ) for m in 1:M )

				Num = 0.5 * ζ^2 * G - ζ * S
				di_Num = di_ζ*ζ * G + 0.5 * ζ^2 * di_G - di_ζ * S - ζ * di_S

				grad[i] = ( di_Num * D - Num * di_D ) / D^2
			end
		end
	end

	s = SobolSeq(N)
	λs = Vector{Float64}[]
	while length(λs) < len_sobol
		λ = next!(s)
		if sum( λ[n] / basis[n].atoms[1].ζ for n in 1:N ) > zero(Float64)
			push!(λs, λ)
		end
	end

	min_couples = pmap(k -> _step_energy_minimization(@inbounds(λs[k]), basis, f, g!, iterations, clean = (rand(1:100) == one(Int64)) ) , 1:len_sobol)
	_, ind = findmin(first.(min_couples))
	cost, minimizer = min_couples[ind]
	return cost, minimizer
end


"""
	energy_errors!( basis::Vector{Mixture{Slater}}, worsts_energy::Vector{Float64}, means_energy::Vector{Float64}, final::Int64, molecules::Vector{Molecule},
					R::Float64, len_sobol::Int64, iterations::Int64, smooth::Float64, save_string::String )
	
	Compute the energy errors until a basis of `final` number of elements.
"""
function energy_errors!( basis::Vector{Mixture{Slater}}, worsts_energy::Vector{Float64}, means_energy::Vector{Float64}, final::Int64, molecules::Vector{Molecule},
						 R::Float64, len_sobol::Int64, iterations::Int64, smooth::Float64, save_string::String )
	@assert length(basis) >= final
	start = length(worsts_energy) + 2
	training_size = length(molecules)
	snapshots = solution.(molecules)

	for n in start:final
		println("Starting online error computations with $n selected elements...")

		wstar = weights_min(basis[1:n])
		mean = 0.0
		worst = -10.0
		best = 10.0
		for k in 1:training_size
			molecule = @inbounds molecules[k]
			cost, λ = energy_minimization(molecule, basis[1:n], wstar, R, len_sobol, iterations, smooth = smooth)
			cost = cost + 0.5 * snapshots[k].atoms[1].ζ^2
			mean += cost
			if cost < best
				best = cost
			end
			if cost > worst
				worst = cost
			end
			println("    minimization $k done with cost $cost.")
		end
		push!(means_energy, mean / training_size)
		push!(worsts_energy, worst)
		save_object(save_string * "energy_errors.jld2", (worsts_energy, means_energy))
		@everywhere GC.gc()
	end
end

