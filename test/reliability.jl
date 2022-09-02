@testset "Reliability" begin
    # set up some irt style test data with correlated responses
    Random.seed!(859345)

    difficulties = randn(10)
    abilities = randn(100)

    m = Matrix{Int}(undef, 100, 10)

    for i in 1:10
        for p in 1:100
            eta = abilities[p] - difficulties[i]
            prob = exp(eta) / (1 + exp(eta))
            m[p, i] = rand(Bernoulli(prob))
        end
    end

    t = ClassicalTestTheory.Test(m)

    @testset "equality of Matrix and Test methods" begin
        @test λ1(m) ≈ λ1(t)
        @test λ2(m) ≈ λ2(t)
        @test λ3(m) ≈ λ3(t)
        @test λ4(m) ≈ λ4(t)
        @test maxλ4(m, method=:bruteforce) ≈ maxλ4(t, method=:bruteforce)
        @test λ5(m) ≈ λ5(t)
        @test λ6(m) ≈ λ6(t)

        @test kr20(m) ≈ kr20(t)
        @test kr21(m) ≈ kr21(t)

        @test glb(m) ≈ glb(t)

    end

    @testset "theoretical guarantees of λ" begin
        @test λ1(t) > 0
        @test λ1(t) < λ3(t) <= λ2(t) <= maxλ4(t, method=:bruteforce)
        @test α(t) == λ3(t)
    end

    # theoretical guarantees of glb
    @testset "theoretical guarantees of glb" begin
        @test maxλ4(t, method=:bruteforce) <= glb(t)
    end
end
