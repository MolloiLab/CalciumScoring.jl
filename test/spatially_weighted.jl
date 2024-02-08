using Test
using CalciumScoring
using Statistics: mean, std

@testset "SpatiallyWeighted" begin
    # 2D
    vol = [
        10 30 135 145 5 150 100 80 9
        9 40 135 130 5 150 190 80 9
        10 100 100 130 5 9 0 80 18
        4 40 135 -10 5 150 100 0 9
        10 0 189 130 5 150 40 80 22
    ]
    calibration = [
        196.569 207.709 177.004 207.147 178.809
        207.359 211.824 0.0 168.913 202.965
        160.564 0.0 208.593 214.862 220.32
        158.685 188.71 226.872 202.076 226.995
        186.996 212.137 185.208 228.175 234.261
    ]
    alg = SpatiallyWeighted()
    answer = 3.7176620417439796
    test = score(vol, calibration, alg)
    @test answer ≈ test
    @test score(vol, mean(calibration), std(calibration), alg) ≈ answer

    # 3D
    vol_ = [
        0 46 58 123 133 104 65 7
        19 25 20 125 163 83 -22 -134
        127 99 65 104 139 57 -51 -102
        123 88 112 140 104 57 46 34
        122 31 71 121 93 83 101 72
        97 24 59 99 68 27 22 52
        0 22 29 72 80 54 28 72
        0 17 -20 15 42 61 80 83
        0 66 14 -3 11 54 110 148
        0 121 127 102 103 93 118 167
    ]
    vol = cat(vol_, vol_, dims=3)
    calibration = [
        153.984 185.388 176.888 168.15 165.752
        155.811 161.226 160.228 183.94 151.235
        149.247 196.509 167.052 158.717 230.236
        174.404 184.419 177.621 149.801 182.128
        161.85 192.528 160.222 178.946 176.578
    ]
    alg = SpatiallyWeighted()
    answer = 1.4580418264134944
    test = score(vol, calibration, alg)
    @test answer ≈ test
    @test score(vol, mean(calibration), std(calibration), alg) ≈ answer
end