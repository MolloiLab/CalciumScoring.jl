abstract type CalciumScore end
struct SpatiallyWeighted <: CalciumScore end

function score(vol, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
	d = Distributions.Normal(μ, σ)

	scaled_array = zeros(size(vol))
	for i in 1:size(vol, 1)
		for j in 1:size(vol, 2)
            if length(size(vol)) == 3
                for z in 1:size(vol, 3)
                    scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
                end
            else
			    scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
            end
		end
	end

    local weighted_arr
    if length(size(vol)) == 3
		kern1 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern2 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern3 = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
		kern = cat(kern1, kern2, kern3, dims=3)
        weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1, 2:end-1]
    else
		kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
	    weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]
    end
    
	return sum(weighted_arr)
end