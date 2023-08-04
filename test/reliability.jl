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

    t = PsychometricTest(m)

    @testset "equality of AbstractMatrix and PsychometricTest methods" begin
        @test λ1(m) ≈ λ1(t)
        @test λ2(m) ≈ λ2(t)
        @test λ3(m) ≈ λ3(t)
        @test λ4(m) ≈ λ4(t)
        @test maxλ4(m, method = :bruteforce) ≈ maxλ4(t, method = :bruteforce)
        @test λ5(m) ≈ λ5(t)
        @test λ6(m) ≈ λ6(t)

        @test kr20(m) ≈ kr20(t)
        @test kr21(m) ≈ kr21(t)

        @test glb(m) ≈ glb(t)
        @test μ(m, 0) ≈ μ(t, 0)
        @test μ(m, 1) ≈ μ(t, 1)
    end

    @testset "λ" begin
        # theoretical guarantees
        @test λ1(t) > 0
        @test λ1(t) < λ3(t) <= λ2(t) <= maxλ4(t, method = :bruteforce)
        @test α(t) == λ3(t)
        @test maxλ4(t, method = :sample, n_samples = 100) <= maxλ4(t, method = :bruteforce)
    end

    @testset "glb" begin
        # theoretical guarantees
        @test maxλ4(t, method = :bruteforce) <= glb(t)
    end

    @testset "μ" begin
        @test_throws ArgumentError μ(m, -1)

        # theoretical guarantees
        for r in 1:10
            @test μ(t, r - 1) <= μ(t, r)
        end

        @test μ(t, 0) ≈ α(t)
        @test μ(t, 1) ≈ λ2(t)
    end
end
