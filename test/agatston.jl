using Revise
using Test
using CalciumScoring
using Unitful: mm, mg

@testset "Agatston" begin
    v1 = ones((4, 2, 2))
    v2 = zeros((4, 2, 2))
    vol = hcat(v1, v2) * 400
    spacing = [0.5, 0.5, 0.5]
    alg = Agatston()

    # Agatston score
    agatston_score, volume_score = score(vol, spacing, alg)
    @test agatston_score ≈ 16 && volume_score == 2

    # Mass score
    mass_cal_factor = 0.00075
    agatston_score, volume_score, mass_score = score(vol, spacing, mass_cal_factor, alg)
    @test agatston_score ≈ 16 && mass_score ≈ 0.6

    # Agatston Unitful
    spacing = [0.5, 0.5, 0.5]mm
    agatston_score, volume_score = score(vol, spacing, alg)
    @test agatston_score ≈ 16 && volume_score == 2mm^3

    # Agatston Mass Unitful
    spacing = [0.5, 0.5, 0.5]mm
    mass_cal_factor = 0.00075mg / mm^3
    agatston_score, volume_score, mass_score = score(vol, spacing, mass_cal_factor, alg)
    @test agatston_score ≈ 16 && mass_score ≈ 0.6mg

    # Various kVs
    hus = [260, 280, 310, 400]
    kVs = [135, 120, 100, 80]
    scores = []
    for i in 1:4
        arr = [
            0 hus[i] hus[i] hus[i] 0 0 0
            0 hus[i] hus[i] hus[i] 0 0 0
            0 hus[i] hus[i] hus[i] 0 0 0
            0 hus[i] hus[i] hus[i] 0 0 0
        ]
        arr = cat(arr, arr, dims=3)
        scr, _ = score(arr, [0.5, 0.5, 0.5], Agatston(); kV=kVs[i])
        push!(scores, scr)
    end
    @test unique(scores) == [12]
end