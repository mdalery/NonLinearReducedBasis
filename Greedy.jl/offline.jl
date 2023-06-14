"""
	projection( mixture::Mixture{Slater}, mixtures::Vector{Mixture{Slater}}, wstar::Dict{CartesianIndex{N}, Float64}, invA::Array{Float64}; clean::Bool = false ) where N

	Compute the distance and weights `λ` of the projection of mixture `mixture` over the set of barycenters between different mixtures `mixtures`.
	Weights ``w^*`` and matrix ``A^{-1}`` needs to be previously computed.
"""
function projection( mixture::Mixture{Slater}, mixtures::Vector{Mixture{Slater}}, wstar::Dict{CartesianIndex{N}, Float64}, invA::Array{Float64}; clean::Bool = false ) where N
	begin
		@assert N isa Int64
		@assert length(mixtures) == N
		@assert size(invA) == (N, N)
	end

	n = length(mixture.weights)
	p = length(wstar)
	coords = collect(keys(wstar))

	c = sum( mixture.weights .* [ atom.r^2 + 2/atom.ζ^2 for atom in mixture.atoms ])

	model = Model()
	@variable(model, w[1:n-1, 1:p-1] >= 0.0)
	@constraint(model, [k = 1:n-1], sum(w[k, :]) <= mixture.weights[k])
	@constraint(model, [kbar = 1:p-1], sum(w[:, kbar]) <= wstar[coords[kbar]])
	@constraint(model, sum(w) >= sum(mixture.weights[1:n-1]) + sum(wstar[c] for c in coords[1:p-1]) - 1.0)

	poly = polyhedron(model, CDDLib.Library())
	vrep_poly = vrep(poly)
	indices = eachindex(points(vrep_poly))
	
	vertex = get(vrep_poly, first(indices))
	w = reshape(vertex, n-1, p-1)
	right = mixture.weights[1:n-1] .- sum(w, dims = 2)
	bottom = [ @inbounds wstar[c] for c in coords[1:p-1] ]' .- sum(w, dims = 1)
	corner = sum(w) - sum(mixture.weights[1:n-1]) - sum( @inbounds wstar[c] for c in coords[1:p-1] ) + 1.0
	w = vcat(hcat(w, right), hcat(bottom, corner))
	bw = -2 .* [ sum( @inbounds w[k, kbar] * (mixture.atoms[k].r * mixtures[i].atoms[coords[kbar][i]].r + 2/(mixture.atoms[k].ζ * mixtures[i].atoms[coords[kbar][i]].ζ)) for k in 1:n, kbar in 1:p ) for i in 1:N ]
	cost = -0.25*bw'*invA*bw + c
	λ = -0.5*invA*bw
	w_min = w
	for index in indices
		vertex = get(vrep_poly, index)

		w = reshape(vertex, n-1, p-1)
		right = mixture.weights[1:n-1] .- sum(w, dims = 2)
		bottom = [ @inbounds wstar[c] for c in coords[1:p-1] ]' .- sum(w, dims = 1)
		corner = sum(w) - sum(mixture.weights[1:n-1]) - sum( @inbounds wstar[c] for c in coords[1:p-1] ) + 1.0
		w = vcat(hcat(w, right), hcat(bottom, corner))
		bw = -2 .* [ sum( @inbounds w[k, kbar] * (mixture.atoms[k].r * mixtures[i].atoms[coords[kbar][i]].r + 2/(mixture.atoms[k].ζ * mixtures[i].atoms[coords[kbar][i]].ζ)) for k in 1:n, kbar in 1:p ) for i in 1:N ]

		cost_w = -0.25*bw'*invA*bw + c
		λ_w = -0.5*invA*bw

		if cost_w < cost
			cost = cost_w
			λ = λ_w
			w_min = w
		end
	end

	if sum( @inbounds λ[i] / mixtures[i].atoms[1].ζ for i in 1:N ) <= 0.0
		error("λ not in Ω_N.")
	end
	
	if clean
		GC.gc()
	end
	return cost, λ, w_min
end


"""
	matrixA( mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64} ) where N

	Compute the matrix ``A``. Weights ``w^*`` needs to be previously computed.
"""
@inline function matrixA( mixtures::Vector{Mixture{Slater}}, w::Dict{CartesianIndex{N}, Float64} ) where N
	begin
		@assert N isa Int64
		@assert length(mixtures) == N
	end

	return [ @inbounds sum( wkbar * ( mixtures[i].atoms[kbar[i]].r * mixtures[j].atoms[kbar[j]].r + 2/( mixtures[i].atoms[kbar[i]].ζ * mixtures[j].atoms[kbar[j]].ζ ) )
				 for (kbar, wkbar) in w ) for i in 1:N, j in 1:N ]
end


"""
	greed!( reduced::Vector{Mixture{Slater}}, N::Int64, snapshots::Vector{Mixture{Slater}}, save_string::String )

	Add the greedy selection of `N` elements in `snapshots` training set to reduced basis `reduced`. Also adds the max and mean errors to `errs_max` and `errs_mean`.
	If `reduced` is not empty, already computed errors are needed.
"""
function greed!( reduced::Vector{Mixture{Slater}}, errs_max::Vector{Float64}, errs_mean::Vector{Float64}, N::Int64, snapshots::Vector{Mixture{Slater}}, save_string::String )
	begin
		@assert N > zero(Int64)
		@assert length(reduced) != 1
		if length(reduced) == 0
			@assert length(errs_max) == length(errs_mean) == 0
		else
			@assert length(errs_max) == length(errs_mean) == length(reduced) - 2
		end
	end
	training_size = length(snapshots)

	start = 1
	if length(reduced) == 0
		h = 0.0
		u1, u2 = @inbounds snapshots[1], snapshots[2]
		for m1 in snapshots
			for m2 in snapshots
				w2 = weights_min([m1, m2])
				h_cur = sq_mixture_distance(m1, m2, w2)
				if h_cur > h
					h = h_cur
					u1, u2 = m1, m2
				end
			end
		end
		push!(reduced, u1, u2)
		start = 3
		println("2 elements selected.")
	end

	for k in start:N
		w = weights_min(reduced)
		A = matrixA(reduced, w)
		errs = pmap(i -> projection(@inbounds(snapshots[i]), reduced, w, inv(A); clean = (rand(1:10) == one(Int64)))[1], 1:training_size)
		err_max, ind_max = findmax(errs)
		err_mean = sum(errs) / training_size
		push!(reduced, @inbounds snapshots[ind_max])
		push!(errs_max, err_max)
		push!(errs_mean, err_mean)

		println("$(length(reduced)) elements selected.")
		save_object(save_string * "/basis.jld2", reduced)
		save_object(save_string * "/projection_errors.jld2", (errs_max, errs_mean))
		println("	saved.")
		
		@everywhere GC.gc()
	end
end
