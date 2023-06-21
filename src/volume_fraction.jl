using Revise
using CalciumScoring

"""
```julia
struct VolumeFraction <: CalciumScore end
```
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the calcium score should be the volume fraction algorithm
"""
struct VolumeFraction <: CalciumScore end


"""
```julia
_percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
```

Find the percentage of calcium within one voxel given the voxel intensity `voxel_intensity`, the Hounsfield unit associated with a known calcium density `hu_calcium`, and the Hounsfield unit associated with background heart tissue `hu_heart_tissue`
"""
function _percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
    return precentage_calcium = (voxel_intensity - hu_heart_tissue) / (hu_calcium - hu_heart_tissue)
end

"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, calculate the number of voxels that are calcified
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            m = vol[i, j]
            percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
            push!(number_calcium_voxels, percent)
        end
    end
    return sum(number_calcium_voxels)
end

"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, and the size of an individual pixel/voxel `voxel_size`, calculate the volume of the calcification
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            m = vol[i, j]
            percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
            push!(number_calcium_voxels, percent)
        end
    end
    return sum(number_calcium_voxels) * voxel_size * density_calcium
end

"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, the size of an individual pixel/voxel `voxel_size`, and the density of the previously mentioned calcium `density_calcium`, calculate the mass of the calcification
"""
function score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            m = vol[i, j]
            percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
            push!(number_calcium_voxels, percent)
        end
    end
    return sum(number_calcium_voxels) * voxel_size * density_calcium
end


"""
```julia
score(vol::AbstractArray, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, calculate the number of voxels that are calcified
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for k in axes(vol, 3)
                m = vol[i, j, k]
                percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
                push!(number_calcium_voxels, percent)
            end
        end
    end
    return sum(number_calcium_voxels)
end

"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, and the size of an individual pixel/voxel `voxel_size`, calculate the volume of the calcification
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for k in axes(vol, 3)
                m = vol[i, j, k]
                percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
                push!(number_calcium_voxels, percent)
            end
        end
    end
    return sum(number_calcium_voxels) * voxel_size * density_calcium
end

"""
```julia
score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
```

Given an input region of interest containing all potential calcium `vol`, the Hounsfield unit associated with a known calcium density `hu_calcium`, the Hounsfield unit associated with background heart tissue `hu_heart_tissue`, the size of an individual pixel/voxel `voxel_size`, and the density of the previously mentioned calcium `density_calcium`, calculate the mass of the calcification
"""
function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
    number_calcium_voxels = []
    for i in axes(vol, 1)
        for j in axes(vol, 2)
            for k in axes(vol, 3)
                m = vol[i, j, k]
                percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
                push!(number_calcium_voxels, percent)
            end
        end
    end
    return sum(number_calcium_voxels) * voxel_size * density_calcium
end
