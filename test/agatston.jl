@testset ExtendedTestSet "Agatston" begin
    @testset ExtendedTestSet "Agatston" begin
        v1 = ones((4, 2, 2))
        v2 = zeros((4, 2, 2))
        vol = hcat(v1, v2) * 400
        spacing = [0.5, 0.5]
        alg = Agatston()
        answer = 16
        @test score(vol, spacing, alg) == answer
    end
end
