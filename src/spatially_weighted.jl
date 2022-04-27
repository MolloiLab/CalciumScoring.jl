# abstract type CalciumScore end
struct SpatiallyWeighted <: CalciumScore end

function score(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)
    μ, σ= mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
		end
	end

	kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]
    
	return sum(weighted_arr)
end

function score(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
		end
	end

	kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]
    
	return sum(weighted_arr)
end

function score(vol::AbstractArray, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			for z in 1:size(vol, 3)
				scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
			end
		end
	end

	weighted_arr = zeros(size(vol))
	for z in 1:size(scaled_array, 3)
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		weighted_arr[:, :, z] = DSP.conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
	end
    
	return sum(weighted_arr)
end

function score(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
			for z in 1:size(vol, 3)
				scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
			end
		end
	end

	weighted_arr = zeros(size(vol))
	for z in 1:size(scaled_array, 3)
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		weighted_arr[:, :, z] = DSP.conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
	end
    
	return sum(weighted_arr)
end