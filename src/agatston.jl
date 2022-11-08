### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 5169922e-50ae-11ed-2236-2f56e5a20360
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise
    using PlutoUI
    using CalciumScoring
end

# ╔═╡ c73efe8b-e72b-48a8-a966-5192d9c99abb
TableOfContents()

# ╔═╡ 3067a00e-be4c-48dc-ac51-1d5e73644f17
"""
```julia
abstract type CalciumScore end
```
Main type for all calcium scoring algorithms in this package
"""
abstract type CalciumScore end

# ╔═╡ 8a560744-3870-4f92-a052-7273259f2bbb
md"""
# Agatston Scoring
"""

# ╔═╡ 40322b74-115f-4785-bbd6-debd834b159c
md"""
## Agatston Score
"""

# ╔═╡ ad02e243-bd61-4a72-913e-1cb901a52cc5
"""
```julia
struct Agatston <: CalciumScore end
```
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the calcium score should be Agatston's algorithm
"""
struct Agatston <: CalciumScore end

# ╔═╡ 728bf5df-c3c7-492e-8735-cb0822e0be20
md"""
## Volume Score
"""

# ╔═╡ ee2a227e-79cf-445a-8c87-45f2738101c2
"""
```julia
function score(vol, spacing, alg::Agatston; kV=120, min_size_mm=1)
```

Given an input `vol` and known pixel/voxel `spacing`, calculate the calcium score via the traditional Agatston scoring technique, as outlined in the [original paper](10.1016/0735-1097(90)90282-T)

Energy (`kV`) specific `threshold`s are determined based on previous [publications](https://doi.org/10.1093/ehjci/jey019)

Returns the Agatston score
"""
function score(vol, spacing, alg::Agatston; kV=120, min_size_mm=1)
    local threshold
    if kV == 80
        threshold = 177
    elseif kV == 100
        threshold = 145
    elseif kV == 120
        threshold = 130
    elseif kV == 135
        threshold = 112
    else
        threshold = Int(round(378 * exp(-0.009 * kV)))
    end
    area_mm = spacing[1] * spacing[2]
    min_size_pixels = Int(round(min_size_mm / area_mm))
    score = 0
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        comp_connect = Int(round(2 * floor(min_size_pixels / 2) + 1))
        if comp_connect > 3
            comp_connect = 3
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        num_non_zero = 0
        number_islands = 0
        slice_score = 0
        num_labels = length(unique(lesion_map))
        for label_idx in 0:num_labels
            if label_idx == 0
                continue
            end

            idxs = findall(x -> x == label_idx, lesion_map)
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

# ╔═╡ fd4ed9d8-c11a-4956-b531-80718047ec90
"""
```julia
function score(vol, spacing, mass_cal_factor, alg::Agatston; kV=120, min_size_mm=1)
```

Given an input `vol`, known pixel/voxel `spacing`, and a mass calibration factor (`mass_cal_factor`) calculate the calcium score via the traditional Agatston scoring technique, as outlined in the [original paper](10.1016/0735-1097(90)90282-T). 

Energy (`kV`) specific `threshold`s are determined based on previous [publications](https://doi.org/10.1093/ehjci/jey019)

Also, converts the Agatston score to a calcium mass score via the `mass_cal_factor`

Returns the Agatston score and the calcium mass score
"""
function score(vol, spacing, mass_cal_factor, alg::Agatston; kV=120, min_size_mm=1)
    local threshold
    if kV == 80
        threshold = 177
    elseif kV == 100
        threshold = 145
    elseif kV == 120
        threshold = 130
    elseif kV == 135
        threshold = 112
    else
        threshold = Int(round(378 * exp(-0.009 * kV)))
    end
    area_mm = spacing[1] * spacing[2]
    slice_thickness = spacing[3]
    min_size_pixels = Int(round(min_size_mm / area_mm))
    mass_score = 0
    score = 0
    for z in axes(vol, 3)
        slice = vol[:, :, z]
        thresholded_slice = slice .> threshold
        max_intensity = maximum(slice)
        if max_intensity < threshold
            continue
        end
        comp_connect = Int(round(2 * floor(min_size_pixels / 2) + 1))
        if comp_connect > 3
            comp_connect = 3
        end
        lesion_map = label_components(thresholded_slice, trues(comp_connect, comp_connect))
        num_non_zero = 0
        number_islands = 0
        mass_slice_score = 0
        slice_score = 0
        num_labels = length(unique(lesion_map))
        for label_idx in 1:(num_labels-1)
            idxs = findall(x -> x == label_idx, lesion_map)
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
            mass_slice_score += mean(intensities)
        end
        plaque_vol = length(findall(x -> x != 0, lesion_map))
        plaque_vol = plaque_vol * area_mm * slice_thickness
        mass_score += mass_slice_score * plaque_vol * mass_cal_factor
        score += slice_score
    end
    return score, mass_score
end

# ╔═╡ Cell order:
# ╠═5169922e-50ae-11ed-2236-2f56e5a20360
# ╠═c73efe8b-e72b-48a8-a966-5192d9c99abb
# ╠═3067a00e-be4c-48dc-ac51-1d5e73644f17
# ╟─8a560744-3870-4f92-a052-7273259f2bbb
# ╟─40322b74-115f-4785-bbd6-debd834b159c
# ╠═ad02e243-bd61-4a72-913e-1cb901a52cc5
# ╟─728bf5df-c3c7-492e-8735-cb0822e0be20
# ╠═ee2a227e-79cf-445a-8c87-45f2738101c2
# ╠═fd4ed9d8-c11a-4956-b531-80718047ec90
