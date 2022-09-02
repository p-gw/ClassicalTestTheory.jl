@testset "Test splitting" begin
    m = zeros(4, 4)
    t = ClassicalTestTheory.Test(m)

    @testset "by UnitRange" begin
        msplit = split(m, 1:2)
        @test length(msplit) == 2
        @test size(first(msplit)) == (4, 2)
        @test size(last(msplit)) == (4, 2)

        tsplit = split(t, 1:2)
        @test length(tsplit) == 2
        @test size(first(tsplit).data) == (4, 2)
        @test size(last(tsplit).data) == (4, 2)
    end

    @testset "by Int" begin
        msplit = split(m, 3)
        @test size(first(msplit)) == (4,)
        @test size(last(msplit)) == (4, 3)

        tsplit = split(t, 3)
        @test size(first(tsplit).data) == (4,)
        @test size(last(tsplit).data) == (4, 3)
    end

    @testset "by Array" begin
        msplit = split(m, [1, 3])
        @test size(first(msplit)) == (4, 2)
        @test size(last(msplit)) == (4, 2)

        tsplit = split(t, [1, 3])
        @test size(first(tsplit).data) == (4, 2)
        @test size(last(tsplit).data) == (4, 2)
    end
end
