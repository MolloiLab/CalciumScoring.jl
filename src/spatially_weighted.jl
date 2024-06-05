using CalciumScoring
using Statistics: mean, std
using Distributions: Normal, cdf
using DSP: conv

"""
## `SpatiallyWeighted`
    
```julia
struct SpatiallyWeighted <: CalciumScore end
```

Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating 
the calcium score should be the Spatially Weighted Calcium Scoring algorithm.
"""
struct SpatiallyWeighted <: CalciumScore end

"""
## `score` (SpatiallyWeighted)

```julia
score(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)

score(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)

score(vol::AbstractArray, calibration, alg::SpatiallyWeighted)

score(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
```

Calculate the calcium score via the Spatially Weighted Calcium Scoring algorithm. This method 
avoids thresholding by weighting each voxel based on a previous calibration (usually 100 mg/cc)
as outlined in the [MESA study](https://doi.org/10.1186/1471-2342-12-14).

The function can be called with different parameters:
1. `vol::AbstractMatrix`, `calibration`, `alg::SpatiallyWeighted` - Applies the algorithm to a 2D volume 
    with a calibration array.
2. `vol::AbstractMatrix`, `μ`, `σ`, `alg::SpatiallyWeighted` - Applies the algorithm to a 2D volume with 
    given mean `μ` and standard deviation `σ`.
3. `vol::AbstractArray`, `calibration`, `alg::SpatiallyWeighted` - Applies the algorithm to a 3D volume with 
    a calibration array.
4. `vol::AbstractArray`, `μ`, `σ`, `alg::SpatiallyWeighted` - Applies the algorithm to a 3D volume with given 
    mean `μ` and standard deviation `σ`.

#### Inputs
- `vol`: input volume containing just the region of interest. This can be a 2D (`AbstractMatrix`) or 
    3D (`AbstractArray`) volume.
- `calibration`: a previous calibration for weighting each voxel.
- `μ`, `σ`: given mean and standard deviation.
- `alg::SpatiallyWeighted`: Spatially Weighted scoring algorithm `SpatiallyWeighted()`.

#### Returns
- `sw_score`: total Spatially Weighted score, providing a continuous and improved measure of atherosclerosis
    compared to existing coronary artery calcium scoring methods.

#### References
[An alternative method for quantifying coronary artery calcification: the multi-ethnic study of atherosclerosis (MESA)](https://doi.org/10.1186/1471-2342-12-14)

[Spatially Weighted Coronary Artery Calcium Score and Coronary Heart Disease Events in the Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1161/CIRCIMAGING.120.011981)
"""
function score(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
    d = Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            scaled_array[i, j] = cdf(d, vol[i, j])
        end
    end

    kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
    weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]
    sw_score = sum(weighted_arr)
    return sw_score
end

function score(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)
    d = Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            scaled_array[i, j] = cdf(d, vol[i, j])
        end
    end

    kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
    weighted_arr = conv(scaled_array, kern)[2:end-1, 2:end-1]
    sw_score = sum(weighted_arr)
    return sw_score
end

function score(vol::AbstractArray, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
    d = Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for z in axes(vol, 3)
                scaled_array[i, j, z] = cdf(d, vol[i, j, z])
            end
        end
    end

    weighted_arr = zeros(size(vol))
    for z in axes(scaled_array, 3)
        kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
        weighted_arr[:, :, z] = conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
    end
    sw_score = sum(weighted_arr)
    return sw_score
end

function score(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
    d = Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for z in axes(vol, 3)
                scaled_array[i, j, z] = cdf(d, vol[i, j, z])
            end
        end
    end

    weighted_arr = zeros(size(vol))
    for z in axes(scaled_array, 3)
        kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
        weighted_arr[:, :, z] = conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
    end
    sw_score = sum(weighted_arr)
    return sw_score
end
