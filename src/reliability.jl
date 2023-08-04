"""
    λ1(m::AbstractMatrix)
    λ1(test::PsychometricTest)
    λ1(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability L₁ described in $GUTTMAN1945
"""
function λ1(m::AbstractMatrix)
    sum_sj = sum(var, eachcol(m))
    st = var(sum(m, dims = 2))
    return 1 - (sum_sj / st)
end

"""
    λ2(scale::AbstractScale)
    λ2(test::PsychometricTest)
    λ2(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability λ₂ described in $GUTTMAN1945
"""
function λ2(m::AbstractMatrix)
    n = size(m, 2)
    C = cov(m)
    zerodiag!(C)
    st = var(sum(m, dims = 2))
    rel = λ1(m) + sqrt(n / (n - 1) * sum(abs2, C)) / st
    return rel
end

"""
    λ3(m::AbstractMatrix)
    λ3(test::PsychometricTest)
    λ3(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability λ₃ described in $GUTTMAN1945
"""
function λ3(m::AbstractMatrix)
    n = size(m, 2)
    return n / (n - 1) * λ1(m)
end

"""
    α(m::AbstractMatrix)
    α(test::PsychometricTest)
    α(test::PsychometricTest, scale::Symbol)

Estimate Cronbach's α. `α` is an alias for [`λ3`](@ref).
"""
const α = λ3

"""
    λ4(m::AbstractMatrix; type::Symbol = :firstlast)
    λ4(test::PsychometricTest; type::Symbol = :firstlast)
    λ4(test::PsychometricTest, scale::Symbol; type::Symbol = :firstlast)

Return the lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

The calculation of λ₄ is based on splitting the test in half.
It is a lower bound of the reliability no matter how the scale is split.

The split of the scale can be controlled by the `type` keyword argument.
The following options are available for `type`:
- `:firstlast`: Split the test by first and last half
- `:oddeven`: Split the test by odd and even indices
- `:random`: Split the test by random indices

To get the maximum lower bound see [`maxλ4`](@ref).
"""
function λ4(m::AbstractMatrix; type::Symbol = :firstlast)
    splits = splithalf(m; type)
    st = var(sum(m, dims = 2))

    s1 = var(sum(splits[1], dims = 2))
    s2 = var(sum(splits[2], dims = 2))

    return 2 * (1 - (s1 + s2) / st)
end

function λ4(m::AbstractMatrix, is)
    splits = split(m, is)
    st = var(sum(m, dims = 2))

    s1 = var(sum(splits[1], dims = 2))
    s2 = var(sum(splits[2], dims = 2))

    return 2 * (1 - (s1 + s2) / st)
end

"""
    maxλ4(m::AbstractMatrix; method = :auto, n_samples = 10_000)
    maxλ4(test::PsychometricTest; method = :auto, n_samples = 10_000)
    maxλ4(test::PsychometricTest, scale::Symbol; method = :auto, n_samples = 10_000)

Calculate the maximum lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

The `method` keyword argument determines the way the bound is estimated. Available options
are:
- `:bruteforce`: Calculate λ₄ for each split-half combination.
- `:sample`: Calculate λ₄ for `n_samples` samples of split-half combinations.
- `:auto` (the default): if the number of items is below 25, `:bruteforce` is applied,
  `:sample` otherwise.

See also [`λ4`](@ref).
"""
function maxλ4(m::AbstractMatrix; method = :auto, n_samples = 10_000)
    if method == :auto
        if size(m, 2) <= 25
            maxλ = _maxλ4_brute_force(m)
        else
            maxλ = _maxλ4_random(m, n_samples)
        end
    elseif method == :bruteforce
        maxλ = _maxλ4_brute_force(m)
    elseif method == :sample
        maxλ = _maxλ4_random(m, n_samples)
    else
        error("Unknown method")
    end

    return maxλ
end

function _maxλ4_brute_force(m::AbstractMatrix)
    n = size(m, 2)
    n_include = ceil(Int, n / 2)
    is = axes(m, 2)

    combs = combinations(is, n_include)
    ncombs = length(combs)

    if ncombs > 1e6
        @info "Brute forcing $(ncombs) combinatinos. Adjust your expectations accordingly..."
    end

    maxλ = maximum(λ4(m, c) for c in combs)

    return maxλ
end

function _maxλ4_random(m::AbstractMatrix, n_samples::Int)
    n = size(m, 2)
    n_include = ceil(Int, n / 2)
    is = axes(m, 2)
    maxλ = maximum(λ4(m, sample(is, n_include, replace = false)) for _ in 1:n_samples)
    return maxλ
end

"""
    λ5(m::AbstractMatrix)
    λ5(test::PsychometricTest)
    λ5(test::PsychometricTest, scale::Symbol)

Return the lower bound estimate of the reliability λ₅ described in $GUTTMAN1945
"""
function λ5(m::AbstractMatrix)
    C = LowerTriangular(cov(m))
    zerodiag!(C)
    sj = sum(abs2, C, dims = 1)
    Cmax = maximum(sj)
    st = var(sum(m, dims = 2))
    return λ1(m) + 2 * sqrt(Cmax) / st
end

"""
    λ6(m::AbstractMatrix)
    λ6(test::PsychometricTest)
    λ6(test::PsychometricTest, scale::Symbol)

Return the lower bound estimate of the reliability λ₆ described in $GUTTMAN1945
"""
function λ6(m::AbstractMatrix)
    C = cov(m)
    Cinv = inv(C)
    smc = 1 .- 1 ./ diag(Cinv)
    return 1 - sum(1 .- smc) / sum(C)
end

"""
    kr20(m::AbstractMatrix)
    kr20(test::PsychometricTest)
    kr20(test::PsychometricTest, scale::Symbol)
"""
function kr20(m::AbstractMatrix)
    n = size(m, 2)
    item_facilities = mean(m, dims = 1)
    item_difficulties = 1 .- item_facilities
    st = var(sum(m, dims = 2))
    return (n / (n - 1)) * ((st - sum(item_facilities .* item_difficulties)) / st)
end

"""
    kr21(m::AbstractMatrix)
    kr21(test::PsychometricTest)
    kr21(test::PsychometricTest, scale::Symbol)
"""
function kr21(m::AbstractMatrix)
    n = size(m, 2)
    item_facilities = mean(m, dims = 1)
    avg_facility = mean(item_facilities)
    avg_difficulty = 1 - avg_facility
    st = var(sum(m, dims = 2))
    return (n / (n - 1)) * ((st - n * avg_difficulty * avg_facility) / st)
end

"""
    glb(m::AbstractMatrix)
    glb(test::PsychometricTest)
    glb(test::PsychometricTest, scale::Symbol)

Return the greatest lower bound estimate (glb) of the reliability as described in
$WOODHOUSE1977
"""
function glb(m::AbstractMatrix)
    n = size(m, 2)

    C = cov(m)
    C̃ = zerodiag(C)
    upr = diag(C)
    lwr = zeros(n)

    model = Model(SCS.Optimizer)

    set_silent(model)
    set_string_names_on_creation(model, false)

    @variable(model, y[1:n])
    @expression(model, A, Symmetric(C̃ + diagm(y)))

    @objective(model, Min, sum(y))
    @constraint(model, lwr .<= y .<= upr)
    @constraint(model, A in PSDCone())

    optimize!(model)

    if termination_status(model) == OPTIMAL
        sum_y = sum(value.(y))
        return (sum(C̃) + sum_y) / sum(C)
    else
        error("something went wrong")
    end
end

"""
    μ(m::AbstractMatrix, r::Int)
    μ(test::PsychometricTest, r::Int)
    μ(test::PsychometricTest, scale::Symbol, r::Int)

Calculate the lower bound of the reliability μ derived in $TENBERGE1978

## Notes
- If `r = 0` then μ is equivalent to Cronbach's alpha.
- If `r = 1` then μ is equivalent to Guttman's λ₂.
"""
function μ(m::AbstractMatrix, r::Int)
    r >= 0 || throw(ArgumentError("r must be non-negative."))

    n = size(m, 2)
    C = BigFloat.(cov(m))
    zerodiag!(C)

    st = var(sum(m, dims = 2))
    p_sum = zero(BigFloat)

    for h in Iterators.reverse(0:r)
        p_h = sum(c -> c^(2^h), C)

        if h == r
            p_h *= n / (n - 1)
        end

        if h == 0
            p_sum = p_sum + p_h
        else
            p_sum = sqrt(p_sum + p_h)
        end
    end

    return Float64(p_sum) / st
end

# generate function definitions for PsychometricTest
for f in [:λ1, :λ2, :λ3, :λ4, :maxλ4, :λ5, :λ6, :kr20, :kr21, :glb, :μ]
    @eval begin
        function $(f)(test::PsychometricTest, args...; kwargs...)
            scales = getscales(test)

            if length(scales) > 0
                rel = [s => $(f)(test, s, args...; kwargs...) for s in keys(scales)]
            else
                rel = $(f)(test, nothing, args...; kwargs...)
            end

            return rel
        end

        function $(f)(test::PsychometricTest, scale::Symbol, args...; kwargs...)
            responses = response_matrix(test, scale)
            return $(f)(responses, args...; kwargs...)
        end

        function $(f)(test::PsychometricTest, ::Nothing, args...; kwargs...)
            responses = response_matrix(test)
            return $(f)(responses, args...; kwargs...)
        end
    end
end
