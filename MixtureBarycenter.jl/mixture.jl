"""
	Mixture{A<:Atom}
	Mixture( atoms::Vector{A}, weights::Vector{Float64} ) where A<:Atom
	
	Object describing a mixture of components `atoms` of type `A` with barycentric weights `weights`.	
"""
struct Mixture{A<:Atom}
	atoms::Vector{A}
	weights::Vector{Float64}
	
	function Mixture( atoms::Vector{A}, weights::Vector{Float64} ) where A<:Atom
		@assert length(atoms) == length(weights)
		@assert are_barycentric(weights)

		new{A}(atoms, weights)
	end
end


"""
	( mixture::Mixture )( x::Float64 )

	The distribution function of the mixture `mixture`.
"""
@inline function ( mixture::Mixture )( x::Float64 )
	sum(mixture.weights .* [ a(x) for a in mixture.atoms ])
end


"""
	weights_min( mixtures::Vector{<:Mixture} )

	Find the barycentric weights ``w^*`` solution to the multi-marginal problem.
	Return a dictionary to mimic a sparse tensor.
"""
function weights_min( mixtures::Vector{<:Mixture} )
	mat = [ sq_atomic_mm(collect(atoms)) for atoms in Iterators.product([ mixture.atoms for mixture in mixtures ]...) ]

	model = Model(GLPK.Optimizer)

	@variable(model, w[ 1:prod(length.(mixture.weights for mixture in mixtures)) ] >= 0.0)
	w = reshape(w, length.(mixture.weights for mixture in mixtures)...)

	N = length(mixtures)
	@constraint(model, [ n=1:N, k=1:length(mixtures[n].weights) ], sum(w[ fill(:, n-1)..., k, fill(:, N-n)... ]) == mixtures[n].weights[k])

	@objective(model, Min, sum(w .* mat))
	optimize!(model)
	w = value.(w)
	Dict( kbar => @inbounds w[kbar] for kbar in CartesianIndices(w) if @inbounds w[kbar] > 0.0 )
end


"""
	sq_mixture_distance( mixture1::M, mixture2::M, w::Dict{CartesianIndex{2}, Float64} ) where M<:Mixture
	
	Compute the squared distance between the mixtures `mixture1` and `mixture2`. Previously computed weights ``w^*`` are required (see weights_min).
"""
@inline function sq_mixture_distance( mixture1::M, mixture2::M, w::Dict{CartesianIndex{2}, Float64} ) where M<:Mixture
	sum( @inbounds wkbar * sq_atomic_distance(mixture1.atoms[kbar[1]], mixture2.atoms[kbar[2]]) for (kbar, wkbar) in w )
end


"""
	sq_mixture_distance( mixture1::M, mixture2::M, w::Dict{CartesianIndex{2}, Float64} ) where M<:Mixture
	
	Compute the distance between the mixtures `mixture1` and `mixture2`. Previously computed weights ``w^*`` are required (see weights_min).
"""
@inline function mixture_distance( mixture1::M, mixture2::M, w::Dict{CartesianIndex{2}, Float64} ) where M<:Mixture
	sqrt(sq_mixture_distance(mixture1, mixture2, w))
end


"""
	barycenter( mixtures::Vector{<:Mixture}, λ::Vector{Float64}, w::Dict{CartesianIndex{N}, Float64} ) where N

	Compute the barycenter between the mixtures `mixtures` with weights `λ`. Previously computed weights ``w^*`` are required (see weights_min).
"""
@inline function barycenter( mixtures::Vector{<:Mixture}, λ::Vector{Float64}, w::Dict{CartesianIndex{N}, Float64} ) where N
	begin
		@assert N isa Int64
		@assert length(mixtures) == length(λ) == N
	end

	atoms = [ barycenter([ @inbounds mixtures[n].atoms[kbar[n]] for n in 1:N ], λ) for kbar in keys(w) ]
	weights = collect(values(w))
	Mixture(atoms, weights)
end
