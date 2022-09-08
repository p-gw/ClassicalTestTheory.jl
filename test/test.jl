@testset "Test" begin
    m = rand(0:1, 10, 10)

    @testset "single scale" begin
        t = ClassicalTestTheory.Test(m)
        s = scales(t)
        @test length(s) == 1
        @test responses(first(s)) == m
        @test id(first(s)) == ""
    end

    @testset "multiple scales" begin
        t = ClassicalTestTheory.Test(
            DichotomousScale("a", m, 1:5),
            DichotomousScale("b", m, 6:10)
        )

        s = scales(t)
        @test length(s) == 2
        @test responses(s[1]) == m[:, 1:5]
        @test responses(s[2]) == m[:, 6:10]

        @test_throws ArgumentError ClassicalTestTheory.Test(
            DichotomousScale("a", m),
            DichotomousScale("a", m)
        )
    end
end
