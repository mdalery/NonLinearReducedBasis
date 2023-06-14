# Offline phase
if N_OFFLINE_SELECT != zero(Int64)
	@info "Performing a basis selection..."
	greed!(basis, worsts_proj, means_proj, N_OFFLINE_SELECT, SNAPSHOTS, SAVE_STRING) # Selection of basis elements
end

# Errors in online phase
if N_MAX_ONLINE < length(basis) && N_MAX_ONLINE != zero(Int64)
	@info "Computing online errors..."
	energy_errors!(basis, worsts_energy, means_energy, N_MAX_ONLINE, MOLECULES, R, LENGTH_SOBOL, ITERATIONS, SMOOTH, SAVE_STRING)
end

# Errors in internal extrapolation phase
if N_MAX_INTER < length(basis) && N_MAX_INTER != zero(Int64)
	@info "Computing internal extrapolation errors..."
	energy_errors!(basis, worsts_energy_inter, means_energy_inter, N_MAX_INTER, MOLECULES_INTER, R, LENGTH_SOBOL, ITERATIONS, SMOOTH, SAVE_STRING * "inter_")
end

# Errors in external extrapolation phase
if N_MAX_EXTER < length(basis) && N_MAX_EXTER != zero(Int64)
	@info "Computing external extrapolation errors..."
	energy_errors!(basis, worsts_energy_exter, means_energy_exter, N_MAX_EXTER, MOLECULES_EXTER, R, LENGTH_SOBOL, ITERATIONS, SMOOTH, SAVE_STRING * "exter_")
end
