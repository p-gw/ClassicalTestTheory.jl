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
        @test lambda1(m) ≈ lambda1(t)
        @test lambda2(m) ≈ lambda2(t)
        @test lambda3(m) ≈ lambda3(t)
        @test lambda4(m) ≈ lambda4(t)
        @test maxlambda4(m, method = :bruteforce) ≈ maxlambda4(t, method = :bruteforce)
        @test lambda5(m) ≈ lambda5(t)
        @test lambda6(m) ≈ lambda6(t)

        @test kr20(m) ≈ kr20(t)
        @test kr21(m) ≈ kr21(t)

        @test glb(m) ≈ glb(t)
        @test mu(m, 0) ≈ mu(t, 0)
        @test mu(m, 1) ≈ mu(t, 1)
    end

    @testset "lambda" begin
        # theoretical guarantees
        @test lambda1(t) > 0
        @test lambda1(t) < lambda3(t) <= lambda2(t) <= maxlambda4(t, method = :bruteforce)
        @test alpha(t) == lambda3(t)
        @test maxlambda4(t, method = :sample, n_samples = 100) <=
              maxlambda4(t, method = :bruteforce)
    end

    @testset "glb" begin
        # theoretical guarantees
        @test maxlambda4(t, method = :bruteforce) <= glb(t)
    end

    @testset "mu" begin
        @test_throws ArgumentError mu(m, -1)

        # theoretical guarantees
        for r in 1:10
            @test mu(t, r - 1) <= mu(t, r)
        end

        @test mu(t, 0) ≈ alpha(t)
        @test mu(t, 1) ≈ lambda2(t)
    end
end
