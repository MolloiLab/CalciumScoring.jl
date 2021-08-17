function agatston_vol(image, spacing, mask, min_vol=nothing, max_vol=nothing)
	binary_mask = mask
	voxel_volume = spacing[1] * spacing[2] * spacing[3]
	
	agatston_score = 0
	calcium_volume = 0
	
	# Find individual lesions (in 3D) so that we can discard too small or too large lesions
	lesion_map = Images.label_components(binary_mask)
	n_lesions = length(unique(lesion_map))
	
	for lesion in 1:(n_lesions - 1)
		lesion_mask = map(x -> x == lesion, lesion_map)
		
		# Ignore too small or too large lesions
		lesion_volume = count(lesion_mask) * voxel_volume
		if ((min_vol != nothing) && (lesion_volume < min_vol))
			continue
		end
		if ((max_vol != nothing) && (lesion_volume > max_vol))
			continue
		end
		
		calcium_volume += lesion_volume
		
		# Calculate Agatston score for this lesion
		slices = sort(unique(lesion_mask)) .+ 1
		for z = 1:slices[end]
			fragment_mask = lesion_mask[z,:,:]
			n_pixels = count(fragment_mask)
			maximum_intensity = maximum(image[z,:,:][fragment_mask])
            if maximum_intensity < 200
                coefficient = 1
			elseif maximum_intensity < 300
				coefficient = 2
			elseif maximum_intensity < 400
				coefficient = 3
			else
				coefficient = 4
			end
			agatston_score += coefficient * n_pixels
		end
	end
	agatston_score *= spacing[1] / 3.0 * spacing[2] * spacing[3]
	return agatston_score, calcium_volume
end
