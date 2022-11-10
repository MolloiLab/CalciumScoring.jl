### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 4d280a50-7a5e-4267-a22e-a372c49658e7
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    using Revise
    using PlutoUI
    using CalciumScoring
end

# ╔═╡ 10173bb4-94dc-49db-b5f9-8206dbfb64ae
TableOfContents()

# ╔═╡ 180d8721-be14-47d2-b50f-6cd60cfac2ed
md"""
# Spatially Weighted Calcium Scoring
"""

# ╔═╡ 91c70333-c03d-4e8d-b382-0b0b20b72884
"""
```julia
struct SpatiallyWeighted <: CalciumScore end
```
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the calcium score should be the Spatially Weighted Calcium Scoring algorithm
"""
struct SpatiallyWeighted <: CalciumScore end

# ╔═╡ 4425fe7f-ebee-4a0e-96eb-cffed0e7f028
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

# ╔═╡ 7d24eba9-9d15-4fcf-90a8-c874e3b2ac35
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

# ╔═╡ 73b7951a-7fa1-493a-a741-312d8c734496
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

# ╔═╡ 65d6d5b5-e73c-4d11-81be-890f02f9740e
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

# ╔═╡ Cell order:
# ╠═4d280a50-7a5e-4267-a22e-a372c49658e7
# ╠═10173bb4-94dc-49db-b5f9-8206dbfb64ae
# ╟─180d8721-be14-47d2-b50f-6cd60cfac2ed
# ╠═91c70333-c03d-4e8d-b382-0b0b20b72884
# ╠═4425fe7f-ebee-4a0e-96eb-cffed0e7f028
# ╠═7d24eba9-9d15-4fcf-90a8-c874e3b2ac35
# ╠═73b7951a-7fa1-493a-a741-312d8c734496
# ╠═65d6d5b5-e73c-4d11-81be-890f02f9740e
