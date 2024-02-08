using Test
using CalciumScoring

@testset "VolumeFraction" begin
    #-- 2D --#
    # Define test data
    vol = [
        400 300
        200 100
        0 0
    ]
    hu_calcium = 400
    hu_heart_tissue = 100

    # Test score function for number of calcified voxels
    s1 = score(vol, hu_calcium, hu_heart_tissue, VolumeFraction())
    @test s1 ≈ 4/3

    # Define voxel_size and density_calcium for further tests
    voxel_size = 0.5
    density_calcium = 1.5

    # Test score function for volume of calcification
    s2 = score(vol, hu_calcium, hu_heart_tissue, voxel_size, VolumeFraction())
    @test s2 ≈ s1 * voxel_size

    # Test score function for mass of calcification
    s3 = score(vol, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, VolumeFraction())
    @test s3 ≈ s1 * voxel_size * density_calcium

    #-- 3D --#
    v = [
        400 300
        200 100
        0 0
    ]
    vol = cat(v, v; dims = 3)

    s1 = score(vol, hu_calcium, hu_heart_tissue, VolumeFraction())
    @test s1 ≈ 8/3

    s2 = score(vol, hu_calcium, hu_heart_tissue, voxel_size, VolumeFraction())
    @test s2 ≈ s1 * voxel_size

    s3 = score(vol, hu_calcium, hu_heart_tissue, voxel_size, density_calcium, VolumeFraction())
    @test s3 ≈ s1 * voxel_size * density_calcium

    zeros1 = score(zeros(3, 3), hu_calcium, 0, VolumeFraction())
    @test zeros1 == 0

    zeros2 = score(zeros(3, 3), hu_calcium, 0, VolumeFraction())
    @test zeros2 == 0
end
