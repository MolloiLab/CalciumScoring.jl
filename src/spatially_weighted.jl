using Revise
using CalciumScoring

"""
```julia
struct SpatiallyWeighted <: CalciumScore end
```
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the calcium score should be the Spatially Weighted Calcium Scoring algorithm
"""
struct SpatiallyWeighted <: CalciumScore end


"""
```julia
score(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)
```

Spatially Weighted Calcium Scoring algorithm. Avoids thresholding by weighting each voxel based on a previous calibration (usually 100 mg/cc)

## Citation
[]
"""
function score(vol::AbstractMatrix, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
    d = Distributions.Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
        end
    end

    kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
    weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]

    return sum(weighted_arr)
end


"""
```julia
score(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)
```

Spatially Weighted Calcium Scoring algorithm. Avoids thresholding by weighting each voxel based on a previous calibration (usually 100 mg/cc)

## Citation
[]
"""
function score(vol::AbstractMatrix, μ, σ, alg::SpatiallyWeighted)
    d = Distributions.Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            scaled_array[i, j] = Distributions.cdf(d, vol[i, j])
        end
    end

    kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
    weighted_arr = DSP.conv(scaled_array, kern)[2:end-1, 2:end-1]

    return sum(weighted_arr)
end


"""
```julia
score(vol::AbstractArray, calibration, alg::SpatiallyWeighted)
```

Spatially Weighted Calcium Scoring algorithm. Avoids thresholding by weighting each voxel based on a previous calibration (usually 100 mg/cc)

## Citation
[]
"""
function score(vol::AbstractArray, calibration, alg::SpatiallyWeighted)
    μ, σ = mean(calibration), std(calibration)
    d = Distributions.Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for z in axes(vol, 3)
                scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
            end
        end
    end

    weighted_arr = zeros(size(vol))
    for z in axes(scaled_array, 3)
        kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
        weighted_arr[:, :, z] = DSP.conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
    end

    return sum(weighted_arr)
end


"""
```julia
score(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
```

Spatially Weighted Calcium Scoring algorithm. Avoids thresholding by weighting each voxel based on a previous calibration (usually 100 mg/cc)

## Citation
[]
"""
function score(vol::AbstractArray, μ, σ, alg::SpatiallyWeighted)
    d = Distributions.Normal(μ, σ)

    scaled_array = zeros(size(vol))
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for z in axes(vol, 3)
                scaled_array[i, j, z] = Distributions.cdf(d, vol[i, j, z])
            end
        end
    end

    weighted_arr = zeros(size(vol))
    for z in axes(scaled_array, 3)
        kern = [0 0.2 0; 0.2 0.2 0.2; 0 0.2 0]
        weighted_arr[:, :, z] = DSP.conv(scaled_array[:, :, z], kern)[2:end-1, 2:end-1]
    end

    return sum(weighted_arr)
end
