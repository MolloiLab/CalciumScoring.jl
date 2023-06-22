using Revise
using CalciumScoring
using Unitful: mm, mg, ustrip, Quantity
using ImageMorphology: label_components

"""
    abstract type CalciumScore end

Main type for all calcium scoring algorithms in this package
"""
abstract type CalciumScore end

"""
    _weight_thresholds(kV, max_intensity)
This code defines a function named `_weight_thresholds` that calculates the weight threshold 
for the Agatston scoring algorithm. The function takes two input parameters: kV and 
max_intensity. The `kV` parameter represents the CT scan's kilovoltage peak (kVp) value, 
and `max_intensity` represents the maximum intensity of a calcified lesion in Hounsfield units (HU).

The function checks the validity of the input `kV` value and ensures that it is one of the accepted 
values: 70, 80, 100, 120, or 135. If the input value does not match any of the accepted values, an 
error message is displayed, and the program is halted.

The function then calculates the weight threshold based on the input `kV` and `max_intensity`. It 
uses a series of if-else statements to determine the weight threshold depending on the kV value and 
the given range of `max_intensity` values.

The function returns the weight threshold as an integer value between 0 and 4.

## Inputs
- `kV`: The kilovoltage peak (kVp) value of the CT scan (Accepted values: 70, 80, 100, 120, or 135)
- `max_intensity`: The maximum intensity of a calcified lesion in Hounsfield units (HU)

## Returns
- `weight`: The weight threshold for the Agatston scoring algorithm (integer value between 0 and 4)

## Reference
- [Gräni, Christoph et al. “Ultra-Low-Dose Coronary Artery Calcium Scoring Using Novel Scoring Thresholds for Low Tube Voltage Protocols—a Pilot Study.” European heart journal cardiovascular imaging 19.12 (2018): 1362–1371. Web.](https://doi.org/10.1093/ehjci/jey019)
"""
function _weight_thresholds(kV, max_intensity)
    @assert (kV == 70 || kV == 80 || kV == 100 || kV == 120 || kV == 135) "kV is $(kV), which is not one of the accepted values \n kV can be either 70, 80, 100, 120 or 135"

    local weight
    if kV == 70
        if max_intensity < 207
            weight = 0
        elseif max_intensity < 318
            weight = 1
        elseif max_intensity < 477
            weight = 2
        elseif max_intensity < 636
            weight = 3
        elseif max_intensity >= 636
            weight = 4
        end
    elseif kV == 80
        if max_intensity < 177
            weight = 0
        elseif max_intensity < 271
            weight = 1
        elseif max_intensity < 408
            weight = 2
        elseif max_intensity < 544
            weight = 3
        elseif max_intensity >= 544
            weight = 4
        end
    elseif kV == 100
        if max_intensity < 145
            weight = 0
        elseif max_intensity < 222
            weight = 1
        elseif max_intensity < 334
            weight = 2
        elseif max_intensity < 446
            weight = 3
        elseif max_intensity >= 446
            weight = 4
        end
    elseif kV == 120
        if max_intensity < 130
            weight = 0
        elseif max_intensity < 199
            weight = 1
        elseif max_intensity < 299
            weight = 2
        elseif max_intensity < 399
            weight = 3
        elseif max_intensity >= 399
            weight = 4
        end
    elseif kV == 135
        if max_intensity < 119
            weight = 0
        elseif max_intensity < 183
            weight = 1
        elseif max_intensity < 274
            weight = 2
        elseif max_intensity < 364
            weight = 3
        elseif max_intensity >= 364
            weight = 4
        end
    end
    return weight
end

"""
    struct Agatston <: CalciumScore end

Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating 
the calcium score should be Agatston's algorithm
"""
struct Agatston <: CalciumScore end

"""
    score(vol, spacing, alg::Agatston; kV=120, min_size_mm=1)
    score(vol, spacing::Array{T}, alg::Agatston; kV=120, min_size_mm=1) where {T<:Unitful.Quantity}
    score(vol, spacing, mass_cal_factor, alg::Agatston; kV=120, min_size_mm=1)
    score(vol, spacing::Array{T}, mass_cal_factor::Unitful.Quantity, alg::Agatston; kV=120, min_size_mm=1) where {T<:Unitful.Quantity}

Calculate the calcium score via the traditional Agatston scoring technique, as outlined in the 
[original paper](10.1016/0735-1097(90)90282-T). Energy (`kV`) specific `threshold`s are determined based on 
previous [publications](https://doi.org/10.1093/ehjci/jey019). Also, it converts the Agatston score to a 
calcium mass score via the `mass_cal_factor` if provided.

## Inputs
- `vol`: input volume containing just the region of interest
- `spacing`: known pixel/voxel spacing (can be an array of `Unitful.Quantity`)
- `alg::Agatston`: Agatston scoring algorithm `Agatston()`
- `mass_cal_factor`: (optional) a factor to calibrate the calculated mass (can be a `Unitful.Quantity`)
- kwargs:
  - `kV=120`: energy of the input CT scan image
  - `min_size=1`: minimum connected component size (see [`label_components`](https://github.com/JuliaImages/Images.jl))

## Returns
- `agatston_score`: total Agatston score
- `volume_score`: total calcium volume via Agatston scoring (can be a `Unitful.Quantity`)
- `abs_mass_score`: (optional) total calcium mass after calibration, returned only when `mass_cal_factor` is provided

## References
[Quantification of coronary artery calcium using ultrafast computed tomography](https://doi.org/10.1016/0735-1097(90)90282-t)

[Ultra-low-dose coronary artery calcium scoring using novel scoring thresholds for low tube voltage protocols—a pilot study ](https://doi.org/10.1093/ehjci/jey019)
"""
function score(vol, spacing, alg::Agatston; kV=120, min_size=1)
    if !(kV in [70, 80, 100, 120, 135])
        error("kV is $(kV), which is not one of the accepted values. kV can be 70, 80, 100, 120, or 135.")
    end
    spacing = ustrip.(spacing)
    threshold = round(378 * exp(-0.009 * kV))
    area = spacing[1] * spacing[2]
    min_size_pixels = div(round(min_size / area), 1)
    comp_connect = Int(min(3, max(1, div(round(2 * div(min_size_pixels, 2) + 1), 1))))
    agatston_score = 0
    volume_score = 0
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        for label_idx in eachindex(unique(lesion_map))
            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end
            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = _weight_thresholds(kV, max_intensity)
            slice_score = num_label_idxs * area * min(weight, 4)
            agatston_score += slice_score
            volume_score += num_label_idxs * spacing[1] * spacing[2] * spacing[3]
        end
    end
    return agatston_score, volume_score
end

function score(vol, spacing::Array{T}, alg::Agatston; kV=120, min_size=1) where {T<:Quantity}
    if !(kV in [70, 80, 100, 120, 135])
        error("kV is $(kV), which is not one of the accepted values. kV can be 70, 80, 100, 120, or 135.")
    end
    spacing = ustrip(spacing)
    threshold = round(378 * exp(-0.009 * kV))
    area = spacing[1] * spacing[2]
    min_size_pixels = div(round(min_size / area), 1)
    comp_connect = Int(min(3, max(1, div(round(2 * div(min_size_pixels, 2) + 1), 1))))
    agatston_score = 0
    volume_score = 0mm^3
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        for label_idx in eachindex(unique(lesion_map))
            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end
            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = _weight_thresholds(kV, max_intensity)
            slice_score = num_label_idxs * area * min(weight, 4)
            agatston_score += slice_score
            volume_score += num_label_idxs * spacing[1]mm * spacing[2]mm * spacing[3]mm
        end
    end
    return agatston_score, volume_score
end


function score(vol, spacing, mass_calibration_factor, alg::Agatston; kV=120, min_size=1)
    if !(kV in [70, 80, 100, 120, 135])
        error("kV is $(kV), which is not one of the accepted values. kV can be 70, 80, 100, 120, or 135.")
    end    
    spacing = ustrip.(spacing)
    threshold = round(378 * exp(-0.009 * kV))
    area = spacing[1] * spacing[2]
    min_size_pixels = div(round(min_size / area), 1)
    comp_connect = Int(min(3, max(1, div(round(2 * div(min_size_pixels, 2) + 1), 1))))
    agatston_score = 0
    volume_score = 0
    attenuations = Float64[]
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        for label_idx in eachindex(unique(lesion_map))
            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end
            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = _weight_thresholds(kV, max_intensity)
            slice_score = num_label_idxs * area * min(weight, 4)
            agatston_score += slice_score
            volume_score += num_label_idxs * spacing[1] * spacing[2] * spacing[3]
            push!(attenuations, intensities...)
        end
    end
    if length(attenuations) == 0
        return agatston_score, volume_score, 0
    else
        rel_mass_score = mean(attenuations) * volume_score
        abs_mass_score = rel_mass_score * mass_calibration_factor
        return agatston_score, volume_score, abs_mass_score
    end
end

function score(vol, spacing::Array{T}, mass_calibration_factor::Quantity, alg::Agatston; kV=120, min_size=1) where {T<:Quantity}
    if !(kV in [70, 80, 100, 120, 135])
        error("kV is $(kV), which is not one of the accepted values. kV can be 70, 80, 100, 120, or 135.")
    end
    spacing = ustrip(spacing)
    threshold = round(378 * exp(-0.009 * kV))
    area = spacing[1] * spacing[2]
    min_size_pixels = div(round(min_size / area), 1)
    comp_connect = Int(min(3, max(1, div(round(2 * div(min_size_pixels, 2) + 1), 1))))
    agatston_score = 0
    volume_score = 0mm^3
    attenuations = []
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        for label_idx in eachindex(unique(lesion_map))
            idxs = findall(x -> x == label_idx, lesion_map)
            num_label_idxs = length(idxs)
            if num_label_idxs < min_size_pixels
                continue
            end
            intensities = slice[idxs]
            max_intensity = maximum(intensities)
            weight = _weight_thresholds(kV, max_intensity)
            slice_score = num_label_idxs * area * min(weight, 4)
            agatston_score += slice_score
            volume_score += num_label_idxs * spacing[1]mm * spacing[2]mm * spacing[3]mm
            push!(attenuations, intensities...)
        end
    end

    local rel_mass_score = 0.0
    if length(attenuations) != 0
        rel_mass_score = mean(attenuations) * volume_score
    end
    abs_mass_score = rel_mass_score * mass_calibration_factor
    return agatston_score, volume_score, abs_mass_score
end
