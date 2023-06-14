"""
	smooth_abs( x::Float64; smooth::Float64 = 1e-8 )

	A smoothed version of the absolute value.
"""
@inline function smooth_abs( x::Float64; smooth::Float64 = 1e-8 )
	if abs(x) < smooth
		return - ((2 * smooth) / π) * cos((π / (2 * smooth)) * x) + smooth
	else
		return abs(x)
	end
end


"""
	smooth_abs( x::Float64; smooth::Float64 = 1e-8 )

	A smoothed version of the sign function.
"""
@inline function smooth_sign( x::Float64; smooth::Float64 = 1e-8 )
	if abs(x) < smooth
		return sin((π / (2 * smooth)) * x)
	else
		return sign(x)
	end
end
