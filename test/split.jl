@testset "Test splitting" begin
    m = zeros(4, 4)

    @testset "by UnitRange" begin
        s = split(m, 1:2)
        @test length(s) == 2
        @test size(first(s)) == (4, 2)
        @test size(last(s)) == (4, 2)
    end

    @testset "by Int" begin
        s = split(m, 3)
        @test size(first(s)) == (4,)
        @test size(last(s)) == (4, 3)
    end

    @testset "by Array" begin
        s = split(m, [1, 3])
        @test size(first(s)) == (4, 2)
        @test size(last(s)) == (4, 2)
    end
end
