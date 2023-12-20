@testset "Reliability" begin
    # set up some irt style test data with correlated responses
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

    @testset "lambda" begin
        # theoretical guarantees
        @test lambda1(m) > 0
        @test lambda1(m) <
              lambda3(m) <=
              lambda2(m) <=
              lambda4(maximum, m, method = :bruteforce).estimate
        @test alpha(m) == lambda3(m)
        @test lambda4(maximum, m, method = :sample, n_samples = 100).estimate <=
              lambda4(maximum, m, method = :bruteforce).estimate
    end

    @testset "glb" begin
        # theoretical guarantees
        @test lambda4(maximum, m, method = :bruteforce).estimate <= glb(m)
    end

    @testset "mu" begin
        @test_throws ArgumentError mu(m, -1)

        # theoretical guarantees
        for r in 1:10
            @test mu(m, r - 1) <= mu(m, r)
        end

        @test mu(m, 0) ≈ alpha(m)
        @test mu(m, 1) ≈ lambda2(m)
    end
end
