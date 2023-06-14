# Examples plot
@info "Generating example plot..."
xs = -4.0:0.001:5.0
solutions = [ solution(Molecule([ 0.8, 1.1 ], [ -r, r ])) for r in [ 0.5, 1.5, 2.5 ] ]

p = plot(xs, [ m.(xs) for m in solutions ],
         linewidth = 5,
		 linestyle = [ :solid :dash :dot ],
         size = (3200, 1000),
         showaxis = :x,
         xticks = -5:1:5,
         tickfontsize = 45,
         legend = :topleft,
         legendfontsize = 45,
		 label = hcat( "solution with r = 0.5",
		 			   "solution with r = 1.5",
		 			   "solution with r = 2.5" ) );

savefig(p, SAVE_STRING_FIGURES * "example_solutions.png")

# Basis plot
@info "Generating basis plot..."
xs = -4.0:0.001:5.0

p = plot(xs, [ m.(xs) for m in basis[1:7] ],
         linewidth = 5,
         size = (3200, 1000),
         showaxis = :x,
         xticks = -5:1:5,
         tickfontsize = 45,
         legend = :topleft,
         legendfontsize = 45,
		 label = hcat( "1st solution selected",
		      		   "2nd solution selected",
					   "3rd solution selected",
					   "4th solution selected",
				 	   "5th solution selected",
					   "6th solution selected",
					   "7th solution selected" ) );

savefig(p, SAVE_STRING_FIGURES * "basis7.png")

# Projection examples
@info "Generating offline projection plots..."
r = 1.266
xs = -4.0:0.001:5.0
molecule = Molecule([ 0.8, 1.1 ], [ -r, r ])
sol = solution(molecule)
Ns = [ 2, 3, 5, 8 ]

for n in Ns
	mixtures = basis[1:n]
	w = weights_min(mixtures)
	A = matrixA(mixtures, w)
	_, λ, _ = projection(sol, mixtures, w, inv(A))
	bar = barycenter(mixtures, λ, w)

	p = plot(xs, [ sol.(xs)  bar.(xs) ],
			 linewidth = 12,
			 linestyle = [ :solid :dash ],
			 size = (3200, 2000),
			 xticks = -5:1:5,
			 tickfontsize = 75,
			 legend = :topleft,
			 label = [ "Solution" "$n selected elements" ],
			 legendfontsize = 75,
			 grid_style = :dot,
			 gridlinewidth = 2.0 );

	   savefig(p, SAVE_STRING_FIGURES * "proj$n.png")
end

# Projection errors plot
@info "Generating offline projection error plot..."
nb_errors = length(worsts_proj) + 1
ymax = log10(worsts_proj[1])
order_max = ceil(Int64, ymax)
ymin = log10(means_proj[end])
order_min = floor(Int64, ymin)
ymax, ymin = 10.0^(ymax + 0.2*(ymax - ymin)), 10.0^(ymin - 0.2*(ymax - ymin))

p = plot(2:nb_errors, [ worsts_proj, means_proj ],
		 linewidth = 10,
		 markershape = [ :star5 :hexagon ],
		 markersize = [40 30],
		 label = [ "Worst MW2 error" "Mean MW2 error" ],
		 legendfontsize = 45,
		 size = (3200, 1800),
		 xticks = 2:nb_errors,
		 tickfontsize = 45,
		 xlabel = "Number of selected elements",
		 xlabelfontsize = 45,
		 grid_style = :dot,
		 gridlinewidth = 2.0,
		 yaxis = :log,
		 yticks = [ 10.0^(order) for order in order_min:order_max ],
		 ylimit = [ ymin, ymax ] );

savefig(p, SAVE_STRING_FIGURES * "projection_errors.png")

# Energy errors plot
@info "Generating energy projection error plot..."
nb_errors = length(worsts_energy) + 1
ymax = log10(worsts_energy[1])
order_max = ceil(Int64, ymax)
ymin = log10(means_energy[end])
order_min = floor(Int64, ymin)
ymax, ymin = 10.0^(ymax + 0.2*(ymax - ymin)), 10.0^(ymin - 0.2*(ymax - ymin))

p = plot(2:nb_errors, [ worsts_energy, means_energy ],
         linewidth = 10,
         markershape = [ :star5 :hexagon ],
         markersize = [40 30],
         label = [ "Worst energy error" "Mean energy error" ],
         legendfontsize = 45,
         size = (3200, 1800),
         xticks = 2:nb_errors,
         tickfontsize = 45,
         xlabel = "Number of selected elements",
         xlabelfontsize = 45,
         grid_style = :dot,
         gridlinewidth = 2.0,
         yaxis = :log,
		 yticks = [ 10.0^(order) for order in order_min:order_max ],
		 ylimit = [ ymin, ymax ] );

savefig(p, SAVE_STRING_FIGURES * "energy_errors.png")

# Internal extrapolation examples
@info "Generating internal energy projection plots..."
r = 0.3
molecule = Molecule([ 0.8, 1.1 ], [ -r, r ])
sol = solution(molecule)
Ns = [ 2, 3, 5, 8 ]
xs = -6.0:0.001:3.0

for n in Ns
	mixtures = basis[1:n]
	w = weights_min(mixtures)
 	_, λ = energy_minimization(molecule, mixtures, w, R, LENGTH_SOBOL, ITERATIONS, smooth = SMOOTH)
 	bar = barycenter(mixtures, λ, w)

 	p = plot(xs, [ sol.(xs)  bar.(xs) ],
			 linewidth = 12,
 			 linestyle = [ :solid :dash ],
 			 size = (3200, 2000),
 			 xticks = -6:1:3,
 			 tickfontsize = 75,
 			 legend = :topleft,
 			 label = [ "Solution" "$n selected elements" ],
 			 legendfontsize = 75,
 			 grid_style = :dot,
 			 gridlinewidth = 2.0 );
 
	savefig(p, SAVE_STRING_FIGURES * "inter$n.png")
end

# Internal extrapolation energy errors plot
@info "Generating internal extrapolation energy projection error plot..."
nb_errors = length(worsts_energy_inter) + 1
ymax = log10(worsts_energy_inter[1])
order_max = ceil(Int64, ymax)
ymin = log10(means_energy_inter[end])
order_min = floor(Int64, ymin)
ymax, ymin = 10.0^(ymax + 0.2*(ymax - ymin)), 10.0^(ymin - 0.2*(ymax - ymin))

p = plot(2:nb_errors, [ worsts_energy_inter, means_energy_inter ],
         linewidth = 10,
         markershape = [ :star5 :hexagon ],
         markersize = [40 30],
         label = [ "Worst energy error" "Mean energy error" ],
         legendfontsize = 45,
         size = (3200, 1800),
         xticks = 2:nb_errors,
         tickfontsize = 45,
         xlabel = "Number of selected elements",
         xlabelfontsize = 45,
         grid_style = :dot,
         gridlinewidth = 2.0,
         yaxis = :log,
		 yticks = [ 10.0^(order) for order in order_min:order_max ],
		 ylimit = [ ymin, ymax ] );

savefig(p, SAVE_STRING_FIGURES * "inter_energy_errors.png")

# External extrapolation examples
@info "Generating external energy projection plots..."
r = 3.5
molecule = Molecule([ 0.8, 1.1 ], [ -r, r ])
sol = solution(molecule)
Ns = [ 2, 3, 5, 8 ]
xs = -3.0:0.001:5.0

for n in Ns
	mixtures = basis[1:n]
	w = weights_min(mixtures)
	_, λ = energy_minimization(molecule, mixtures, w, R, LENGTH_SOBOL, ITERATIONS, smooth = SMOOTH)
	bar = barycenter(mixtures, λ, w)

	p = plot(xs, [ sol.(xs)  bar.(xs) ],
			 linewidth = 12,
			 linestyle = [ :solid :dash ],
			 size = (3200, 2000),
			 xticks = -3:1:5,
			 tickfontsize = 75,
			 legend = :topleft,
			 label = [ "Solution" "$n selected elements" ],
			 legendfontsize = 75,
			 grid_style = :dot,
			 gridlinewidth = 2.0 );

	   savefig(p, SAVE_STRING_FIGURES * "exter$n.png")
end

# External extrapolation energy errors plot
@info "Generating external extrapolation energy projection error plot..."
nb_errors = length(worsts_energy_exter) + 1
ymax = log10(worsts_energy_exter[1])
order_max = ceil(Int64, ymax)
ymin = log10(means_energy_exter[end])
order_min = floor(Int64, ymin)
ymax, ymin = 10.0^(ymax + 0.2*(ymax - ymin)), 10.0^(ymin - 0.2*(ymax - ymin))

p = plot(2:nb_errors, [ worsts_energy_exter, means_energy_exter ],
         linewidth = 10,
         markershape = [ :star5 :hexagon ],
         markersize = [40 30],
         label = [ "Worst energy error" "Mean energy error" ],
         legendfontsize = 45,
         size = (3200, 1800),
         xticks = 2:nb_errors,
         tickfontsize = 45,
         xlabel = "Number of selected elements",
         xlabelfontsize = 45,
         grid_style = :dot,
         gridlinewidth = 2.0,
         yaxis = :log,
		 yticks = [ 10.0^(order) for order in order_min:order_max ],
		 ylimit = [ ymin, ymax ] );

savefig(p, SAVE_STRING_FIGURES * "exter_energy_errors.png")

# Energy profile
@info "Generating energy heatmap plot..."
r = 2.15
molecule = Molecule([ 0.8, 1.1 ], [ -r, r ])
xs = -2.0:0.01:4.0
ys = -4.0:0.01:2.0
mixtures = basis[1:2]
w = weights_min(mixtures)
zs = [ try; energy([x, y], molecule, mixtures, w); catch; NaN; end for x in xs, y in ys ];

p = heatmap(xs, ys, zs',
			size = (3200, 1800),
			tickfontsize = 45,
			colorbar_ticks = -0.5:0.25:1.0,
			colorbar_tickfontsize = 45,
			fillcolor = cgrad(:matter, alpha = 1.0),
			clims = (-0.6, 1.0) );


savefig(p, SAVE_STRING_FIGURES * "energy_heatmap.png")
