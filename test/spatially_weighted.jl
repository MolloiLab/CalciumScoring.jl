@testset ExtendedTestSet "SpatiallyWeighted" begin
    @testset ExtendedTestSet "SpatiallyWeighted 2D" begin
        vol = [
            10 30 135 145 5 150 100 80 9
            9 40 135 130 5 150 190 80 9
            10 100 100 130 5 9 0 80 18
            4 40 135 -10 5 150 100 0 9
            10 0 189 130 5 150 40 80 22
        ]
        calibration = [
            196.569  207.709  177.004  207.147  178.809
            207.359  211.824    0.0    168.913  202.965
            160.564    0.0    208.593  214.862  220.32
            158.685  188.71   226.872  202.076  226.995
            186.996  212.137  185.208  228.175  234.261
        ]
        alg = SpatiallyWeighted()
        answer = 3.7176620417439796
        test = score(vol, calibration, alg)
        @test answer ≈ test
    end
    @testset ExtendedTestSet "SpatiallyWeighted 3D" begin
        vol1 = [
            10 30 135 145 5 150 100 80 9
            9 40 135 130 5 150 190 80 9
            10 100 100 130 5 9 0 80 18
            4 40 135 -10 5 150 100 0 9
            10 0 189 130 5 150 40 80 22
        ]
        vol = cat(vol1, vol1, dims=3)
        calibration1 = [
            196.569  207.709  177.004  207.147  178.809
            207.359  211.824    0.0    168.913  202.965
            160.564    0.0    208.593  214.862  220.32
            158.685  188.71   226.872  202.076  226.995
            186.996  212.137  185.208  228.175  234.261
        ]
        calibration = cat(calibration1, calibration1, dims=3)
        alg = SpatiallyWeighted()
        answer = 14.71253963876056
        test = score(vol, calibration, alg)
        @test answer ≈ test
    end
end