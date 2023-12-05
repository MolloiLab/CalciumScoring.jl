using Random

function create_calcium_calibration_phantom(
    phantom_size::Tuple{Int64, Int64}, 
	calcium_densities;
    noise_mean = 0, noise_std = 50, energy = 120
)

    # Adjust calcium intensities based on energy
    function adjust_for_energy(intensity, energy)
        if energy == 135
            return intensity * 0.82  # Adjust these factors as needed
        elseif energy == 80
            return intensity * 1.14  # Adjust these factors as needed
        else
            return intensity
        end
    end

    phantom = fill(-1024.0, phantom_size...)
    masks = zeros(phantom_size..., length(calcium_densities))  # Initialize a 3D mask array with zeros

    center_x, center_y = div(phantom_size[1], 2), div(phantom_size[2], 2)
    soft_tissue_radius = div(256, 2)

    # Function to check if a point is inside a circle
    is_inside_circle(x, y, center_x, center_y, radius) = (x - center_x)^2 + (y - center_y)^2 <= radius^2

    # Draw the soft tissue circle
    for x in 1:phantom_size[1]
        for y in 1:phantom_size[2]
            if is_inside_circle(x, y, center_x, center_y, soft_tissue_radius)
                phantom[x, y] = 0
            end
        end
    end

    # Calcium insert details
    calcium_intensities = 1.4 .* 1e3 .* adjust_for_energy.(calcium_densities, energy)
    calcium_insert_radius = div(25, 2)
    angle_increment = 2 * π / length(calcium_intensities)
    distance_to_center = 0.6 * (soft_tissue_radius - calcium_insert_radius * 2)

    # Place the adjusted calcium inserts and create individual masks
    for (idx, intensity) in enumerate(calcium_intensities)
        angle = idx * angle_increment
        insert_center_x = center_x + distance_to_center * cos(angle)
        insert_center_y = center_y + distance_to_center * sin(angle)

        for x in 1:phantom_size[1]
            for y in 1:phantom_size[2]
                if is_inside_circle(x, y, insert_center_x, insert_center_y, calcium_insert_radius)
                    phantom[x, y] = intensity
                    masks[x, y, idx] = 1
                end
            end
        end
    end

    # Add Gaussian noise to the phantom
    noise = noise_mean .+ noise_std .* randn(phantom_size...)
    phantom .+= noise

    return phantom, masks
end

function create_calcium_calibration_phantom(
    phantom_size::Tuple{Int64, Int64, Int64},
	calcium_densities; 
    noise_mean = 0, noise_std = 50, energy = 120
)
    phantom_slices = []
    mask_slices = []

    for _ in 1:phantom_size[3]
        p, m = create_calcium_calibration_phantom(
			phantom_size[1:2], calcium_densities;
			noise_mean=noise_mean, noise_std=noise_std, energy=energy)
        push!(phantom_slices, p)
        push!(mask_slices, m)
    end

    phantom_3d = cat(phantom_slices..., dims=3)
    masks_3d = cat(mask_slices..., dims=4)
    masks_3d = permutedims(masks_3d, (1, 2, 4, 3))

    return phantom_3d, masks_3d
end

function create_calcium_measurement_phantom(
    phantom_size::Tuple{Int64, Int64},
    calcium_densities;
    noise_mean = 0, noise_std = 50, energy = 120
)

    # Adjust calcium intensities based on energy
    function adjust_for_energy(intensity, energy)
        if energy == 135
            return intensity * 0.82
        elseif energy == 80
            return intensity * 1.14
        else
            return intensity
        end
    end

    phantom = fill(-1024.0, phantom_size...)
    masks = zeros(phantom_size..., length(calcium_densities))

    center_x, center_y = div(phantom_size[1], 2), div(phantom_size[2], 2)
    soft_tissue_radius = div(256, 2)

    # Function to check if a point is inside a circle
    is_inside_circle(x, y, center_x, center_y, radius) = (x - center_x)^2 + (y - center_y)^2 <= radius^2

    # Draw the soft tissue circle
    for x in 1:phantom_size[1]
        for y in 1:phantom_size[2]
            if is_inside_circle(x, y, center_x, center_y, soft_tissue_radius)
                phantom[x, y] = 0
            end
        end
    end

    # Adjusted calcium insert details
    calcium_intensities = 1e3 .* adjust_for_energy.(calcium_densities, energy)
    calcium_insert_radius = div(15, 2)  # Reduced size for the measurement phantom inserts
    
	# Adjusted calcium insert details
	num_inserts = length(calcium_densities)
	angle_increment = 2 * π / num_inserts
	distance_to_center = 50  # Adjusted distance to center for measurement inserts
	
	insert_positions = [
	    (center_x + distance_to_center * cos(idx * angle_increment),
	     center_y + distance_to_center * sin(idx * angle_increment))
	    for idx in 1:num_inserts
	]

    # Place the calcium inserts and create individual masks
    for (idx, (intensity, (insert_center_x, insert_center_y))) in enumerate(zip(calcium_intensities, insert_positions))
        for x in 1:phantom_size[1]
            for y in 1:phantom_size[2]
                if is_inside_circle(x, y, insert_center_x, insert_center_y, calcium_insert_radius)
                    phantom[x, y] = intensity
                    masks[x, y, idx] = 1
                end
            end
        end
    end

    # Add Gaussian noise to the phantom
    noise = noise_mean .+ noise_std .* randn(phantom_size...)
    phantom .+= noise

    return phantom, masks
end

function create_calcium_measurement_phantom(
    phantom_size::Tuple{Int64, Int64, Int64},
    calcium_densities; 
    noise_mean = 0, noise_std = 50, energy = 120
)
    phantom_slices = []
    mask_slices = []

    for _ in 1:phantom_size[3]
        p, m = create_calcium_measurement_phantom(
            phantom_size[1:2], calcium_densities;
            noise_mean=noise_mean, noise_std=noise_std, energy=energy)
        push!(phantom_slices, p)
        push!(mask_slices, m)
    end

    phantom_3d = cat(phantom_slices..., dims=3)
    masks_3d = cat(mask_slices..., dims=4)
    masks_3d = permutedims(masks_3d, (1, 2, 4, 3))

    return phantom_3d, masks_3d
end