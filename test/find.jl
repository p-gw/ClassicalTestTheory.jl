@testset "find" begin
    # setup data
    n_items = 4
    n_persons = 100
    ρ = 0.8

    μ = zeros(n_items)
    Σ = fill(ρ, n_items, n_items)

    for i in 1:n_items
        Σ[i, i] = 1
    end

    dist = MvNormal(μ, Σ)

    m = rand(dist, n_persons)' .< 0
    m_extended = hcat(m, zeros(n_persons))

    # test
    @test_throws ArgumentError find(m, n_items + 2)
    @test size(find(m, 2, progress = false)) == (n_persons, 2)
    @test size(find(m, 1, progress = false)) == (n_persons, 1)

    @test find(m_extended, n_items, progress = false) == m
    @test find(m, 2, progress = true) == find(m, 2, progress = false)
end

