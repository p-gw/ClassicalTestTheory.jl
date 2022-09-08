@testset "Scales" begin
    @testset "DichotomousScale" begin
        @testset "Matrix{Bool}" begin
            m = rand([true, false], 10, 10)
            scale = DichotomousScale("a", m)
            @test id(scale) == "a"
            @test responses(scale) == m
            @test scores(scale) == vec(sum(m, dims=2))
            @test cov(scale) == cov(m)

            scale = DichotomousScale("a", m, 1:5)
            @test responses(scale) == m[:, 1:5]

            scale = DichotomousScale("a", m, 6)
            @test responses(scale) == reshape(m[:, 6], :, 1)
        end

        @testset "Matrix{Int}" begin
            m = rand(0:1, 10, 10)
            scale = DichotomousScale("b", m)
            @test id(scale) == "b"
            @test responses(scale) == m
            @test scores(scale) == vec(sum(m, dims=2))
            @test cov(scale) == cov(m)

            scale = DichotomousScale("b", m, 1:4)
            @test responses(scale) == m[:, 1:4]

            scale = DichotomousScale("b", m, 10)
            @test responses(scale) == reshape(m[:, 10], :, 1)

            m = rand(1:5, 10, 10)
            @test_throws InexactError DichotomousScale("", m)
        end

        # TODO: DataFrame
        @testset "DataFrame" begin end
    end
end
