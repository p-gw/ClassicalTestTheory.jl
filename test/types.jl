@testset "Types" begin
    @testset "Matrix methods" begin
        m = [1 0; 0 1; 1 1]

        @test scores(m) == [1, 1, 2]

        @test responses(m) == m
        @test responses(m, 1) == [1, 0, 1]
        @test responses(m, 2) == [0, 1, 1]
        @test_throws BoundsError responses(m, 3)

        @test nitems(m) == 2
        @test npersons(m) == 3
    end

    @testset "Test" begin
        m = [1 0; 0 1; 1 1]
        t = ClassicalTestTheory.Test(m)

        @test scores(t) == [1, 1, 2]
        @test responses(t) == m
        @test responses(t, 1) == [1, 0, 1]
        @test_throws BoundsError responses(t, 5)

        @test cov(t) == t.itemcov == cov(m)
    end

    @testset "SubTest" begin
        m = [1 0; 0 1; 1 1]
        t = ClassicalTestTheory.Test(m)

        @testset "is::Int" begin
            st = SubTest(t, 1)
            @test scores(st) == [1, 0, 1]
            @test responses(st) == [1, 0, 1]
            @test responses(st, 1) == [1, 0, 1]
            @test_throws BoundsError responses(st, 2)
        end

        @testset "is::UnitRange" begin
            st = SubTest(t, 1:2)
            @test scores(st) == [1, 1, 2]
            @test responses(st) == m
            @test responses(st, 1) == [1, 0, 1]
            @test_throws BoundsError responses(st, 3)
        end

        @testset "is::Vector" begin
            st = SubTest(t, [1, 2])
            @test scores(st) == [1, 1, 2]
            @test responses(st) == m
            @test responses(st, 1) == [1, 0, 1]
            @test_throws BoundsError responses(st, 3)
        end
    end
end
