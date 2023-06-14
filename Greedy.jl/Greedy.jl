############################################################################################
#####                                                                                  ##### 
#####                                                                                  ##### 
#####     Greedy.jl                                                                    #####
#####                                                                                  #####
#####                                                                                  #####
#####     Julia module for mixture barycenter projection and energy minimization.      #####
#####                                                                                  #####
#####                                                                                  ##### 
#####     Made with julia 1.9.0                                                        #####
#####                                                                                  ##### 
#####                                                                                  ##### 
############################################################################################

module Greedy

using Distributed,
	  JLD2,
	  LinearAlgebra,
	  JuMP,
	  Polyhedra,
	  CDDLib,
	  Sobol,
	  Optim,
	  Main.MixtureBarycenter

include("util.jl")

include("solutions.jl")
export Molecule,
	   solution

include("offline.jl")
export matrixA,
	   projection,
	   greed!

include("online.jl")
export energy,
	   senergy,
	   di_energy,
	   di_senergy,
	   energy_minimization,
	   energy_errors!

end
