abstract type CalciumScore end

# struct Agatston{T1<:AbstractArray,T2<:AbstractFloat} <: CalciumScore
#     vol::T1
#     spacing::T2
# end

# function Agatston(
#     vol::AbstractArray,
#     spacing::AbstractArray
# )
#     return Agatston(vol, spacing)
# end

struct Agatston <: CalciumScore end

function score(vol, spacing, alg::Agatston; threshold=130, min_size_mm=1)
	area_mm = spacing[1] * spacing[2]
	min_size_pixels = Int(round(min_size_mm / area_mm))
	num_slices = size(vol, 3)
	score = 0
	for z in range(1, num_slices)
		slice = vol[:, :, z]
		thresholded_slice = slice .> threshold
		max_intensity = maximum(slice)
		if max_intensity < threshold
			continue
		end
		comp_connect = trues(min_size_pixels + 1, min_size_pixels + 1)
		lesion_map = label_components(thresholded_slice, comp_connect)
		num_non_zero = 0
		number_islands = 0
		slice_score = 0
		num_labels = length(unique(lesion_map))
		for label_idx in 0:num_labels
			if label_idx == 0
				continue
			end

			idxs = findall(x->x==label_idx, lesion_map)
			num_label_idxs = length(idxs)
			if num_label_idxs < min_size_pixels
				continue
			end

			intensities = slice[idxs]
			max_intensity = maximum(intensities)
			weight = floor(max_intensity / 100)
			if weight > 4
				weight = 4.0
			end
			num_non_zero_in_island = num_label_idxs
			slice_score += num_non_zero_in_island * area_mm * weight
			num_non_zero += num_non_zero_in_island
			number_islands += 1
		end
		score += slice_score
	end
	return score
end