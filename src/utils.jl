"""
	mask_elements(img, threshold, connectivity)
Remove every element that is below the `threshold` and
adjacent to more than the `connectivity` amount of other 
elements similarly below the threshold. Returns a mask where 
all elements are 1 for elements that should be kept based on
`threshold` and `connectivity` criteria.

# Arguments
- img: array which will be used to determine mask
- threshold: cutoff value for elements that should be 
    removed
- connectivity: minimum number of adjacent elements that are 
    required to remove an element
"""
function mask_elements(img, threshold, connectivity)
    mask = img .< threshold
    label = Images.label_components(mask)
    n = length(unique(label))
    indices = []
    for lbl in 0:(n - 1)
        lbl_mask = map(x -> x == lbl, label)
        num = count(lbl_mask)
        if num < connectivity
            rmv_indices = findall(x -> x == lbl, label)
            push!(indices, rmv_indices)
        end
    end
    for i in 1:length(indices)
        idx = indices[i]
        mask[idx] .= 0
    end
    return iszero.(mask)
end

"""

"""
function mass_calibration(
    dcm_array, center_large_LD, center, cal_rod_slice, rows, cols, spacing
)
    center_LD = center_large_LD
    dist_x = abs(center_LD[1] - center[1])
    dist_y = abs(center_LD[2] - center[2])

    if dist_x == 0
        mass_center_x = center[1]
        if center_LD[2] > center[2]
            mass_center_y = round(center[1] - round(23 / spacing[1], RoundUp), RoundUp)
        else
            mass_center_y = round(center[1] + round(23 / spacing[1], RoundUp), RoundUp)
        end
    elseif dist_y == 0
        mass_center_y = center[2]
        if center_LD[1] > center[1]
            mass_center_x = round(center[1] - round(23 / spacing[1], RoundUp), RoundUo)
        else
            mass_center_x = round(center[0] + round(23 / spacing[1], RoundUp), RoundUp)
        end

    else
        mass_angle = atan(dist_y / dist_x)
        dist_x = (23 / spacing[1]) * cos(mass_angle)
        dist_y = (23 / spacing[1]) * sin(mass_angle)

        if (center_LD[1] < center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] < center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] + dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] < center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] + dist_y, RoundUp)
        elseif (center_LD[1] > center[1] && center_LD[2] > center[2])
            mass_center_x = round(center[1] - dist_x, RoundUp)
            mass_center_y = round(center[2] - dist_y, RoundUp)
        end
    end

    mass_cal_center = [mass_center_y, mass_center_x]
    x_distance = abs(center[1] - mass_cal_center[2])
    angled_distance = sqrt(
        (center[1] - mass_cal_center[2])^2 + (center[2] - mass_cal_center[2])^2
    )
    angle_0_200HA = acos(x_distance / angled_distance) * 180 / π
    mask_0HU = Phantoms.create_circular_mask(
        cols, rows, mass_cal_center, Int(round(6.9 / spacing[1]))
    )
    masked_0HU = mask_0HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count = length(findall(x -> x != 0, masked_0HU))
    mean_0HU = sum(masked_0HU) / nonzero_count
    std_0HU = 0

    for voxel in vec(masked_0HU)
        if voxel != 0
            std_0HU += (voxel - mean_0HU)^2
        end
    end

    std_0HU = sqrt(std_0HU / nonzero_count)
    mask_200HU = Phantoms.create_circular_mask(
        cols, rows, (center[2], center[1]), Int(round(6.9 / spacing[1]))
    )
    masked_200HU = mask_200HU .* dcm_array[:, :, cal_rod_slice]
    nonzero_count_200HU = length(findall(x -> x != 0, masked_200HU))
    mean_200HU = sum(masked_200HU) / nonzero_count_200HU
    mass_cal_factor = 0.2 / (mean_200HU - mean_0HU)
    water_rod_metrics = mean_0HU, std_0HU

    return mass_cal_factor, angle_0_200HA, water_rod_metrics
end
