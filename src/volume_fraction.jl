using Revise
using CalciumScoring

"""
    struct VolumeFraction <: CalciumScore end
Lets Julia know (via multiple dispatch) that the algorithm of choice when calculating the 
calcium score should be the volume fraction algorithm
"""
struct VolumeFraction <: CalciumScore end


"""
    _percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)

Find the percentage of calcium within one voxel given the voxel intensity `voxel_intensity`, 
the Hounsfield unit associated with a known calcium density `hu_calcium`, and the Hounsfield 
unit associated with background heart tissue `hu_heart_tissue`
"""
function _percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
    return precentage_calcium = (voxel_intensity - hu_heart_tissue) / (hu_calcium - hu_heart_tissue)
end

"""
    score(vol::AbstractArray, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
    score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
    score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)
    score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, alg::VolumeFraction)
    score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, alg::VolumeFraction)
    score(vol::AbstractMatrix, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, alg::VolumeFraction)

Calculate the number of voxels that are calcified, the volume of the calcification, or the mass of the calcification 
depending on the provided arguments.

#### Inputs
- `vol`: An input region of interest containing all potential calcium.
- `hu_calcium`: The Hounsfield unit associated with a known calcium density.
- `hu_heart_tissue`: The Hounsfield unit associated with background heart tissue.
- `voxel_size`: (optional) The size of an individual pixel/voxel.
- `density_calcium`: (optional) The density of the previously mentioned calcium.
- `alg::VolumeFraction`: The algorithm to be used for the calculation.

#### Returns
- When called with `vol`, `hu_calcium`, `hu_heart_tissue`, and `alg`, this function calculates and returns 
    the number of voxels that are calcified.
- When called with `vol`, `hu_calcium`, `hu_heart_tissue`, `voxel_size`, and `alg`, this function calculates 
    and returns the volume of the calcification.
- When called with `vol`, `hu_calcium`, `hu_heart_tissue`, `voxel_size`, `density_calcium`, and `alg`, this 
    function calculates and returns the mass of the calcification.
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
