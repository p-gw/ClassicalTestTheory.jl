"""
    lambda1(m::AbstractMatrix)
    lambda1(test::PsychometricTest)
    lambda1(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability L₁ described in $GUTTMAN1945
"""
function lambda1(m::AbstractMatrix)
    sum_sj = sum(var, eachcol(m))
    st = var(sum(m, dims = 2))
    return 1 - (sum_sj / st)
end

struct L1 <: ReliabilityMeasure end
(method::L1)(data) = lambda1(data)

"""
    lambda2(scale::AbstractScale)
    lambda2(test::PsychometricTest)
    lambda2(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability lambda₂ described in $GUTTMAN1945
"""
function lambda2(m::AbstractMatrix)
    n = size(m, 2)
    C = cov(m)
    zerodiag!(C)
    st = var(sum(m, dims = 2))
    rel = lambda1(m) + sqrt(n / (n - 1) * sum(abs2, C)) / st
    return rel
end

struct L2 <: ReliabilityMeasure end
(method::L2)(data) = lambda2(data)

"""
    lambda3(m::AbstractMatrix)
    lambda3(test::PsychometricTest)
    lambda3(test::PsychometricTest, scale::Symbol)

Calculate the lower bound estimate of the reliability lambda₃ described in $GUTTMAN1945
"""
function lambda3(m::AbstractMatrix)
    n = size(m, 2)
    return n / (n - 1) * lambda1(m)
end

struct L3 <: ReliabilityMeasure end
(method::L3)(data) = lambda3(data)

"""
    alpha(m::AbstractMatrix)
    alpha(test::PsychometricTest)
    alpha(test::PsychometricTest, scale::Symbol)

Estimate Cronbach's alpha. `alpha` is an alias for [`lambda3`](@ref).
"""
const alpha = lambda3

struct Alpha <: ReliabilityMeasure end
(method::Alpha)(data) = alpha(data)

"""
    lambda4(m::AbstractMatrix; type::Symbol = :firstlast)
    lambda4(test::PsychometricTest; type::Symbol = :firstlast)
    lambda4(test::PsychometricTest, scale::Symbol; type::Symbol = :firstlast)

Return the lower bound estimate of the reliability lambda₄ described in $GUTTMAN1945

The calculation of lambda₄ is based on splitting the test in half.
It is a lower bound of the reliability no matter how the scale is split.

The split of the scale can be controlled by the `type` keyword argument.
The following options are available for `type`:
- `:firstlast`: Split the test by first and last half
- `:oddeven`: Split the test by odd and even indices
- `:random`: Split the test by random indices

To get the maximum lower bound see [`maxlambda4`](@ref).
"""
function lambda4(m::AbstractMatrix; type::Symbol = :firstlast)
    splits = splithalf(m; type)
    st = var(sum(m, dims = 2))

    s1 = var(sum(splits[1], dims = 2))
    s2 = var(sum(splits[2], dims = 2))

    return 2 * (1 - (s1 + s2) / st)
end

function lambda4(m::AbstractMatrix, is)
    splits = split(m, is)
    st = var(sum(m, dims = 2))

    s1 = var(sum(splits[1], dims = 2))
    s2 = var(sum(splits[2], dims = 2))

    return 2 * (1 - (s1 + s2) / st)
end


"""
    maxlambda4(m::AbstractMatrix; method = :auto, n_samples = 10_000)
    maxlambda4(test::PsychometricTest; method = :auto, n_samples = 10_000)
    maxlambda4(test::PsychometricTest, scale::Symbol; method = :auto, n_samples = 10_000)

Calculate the maximum lower bound estimate of the reliability lambda₄ described in $GUTTMAN1945

The `method` keyword argument determines the way the bound is estimated. Available options
are:
- `:bruteforce`: Calculate lambda₄ for each split-half combination.
- `:sample`: Calculate lambda₄ for `n_samples` samples of split-half combinations.
- `:auto` (the default): if the number of items is below 25, `:bruteforce` is applied,
  `:sample` otherwise.

See also [`lambda4`](@ref).
"""
function maxlambda4(m::AbstractMatrix; method = :auto, n_samples = 10_000)
    if method == :auto
        if size(m, 2) <= 25
            maxlambda = _maxlambda4_brute_force(m)
        else
            maxlambda = _maxlambda4_random(m, n_samples)
        end
    elseif method == :bruteforce
        maxlambda = _maxlambda4_brute_force(m)
    elseif method == :sample
        maxlambda = _maxlambda4_random(m, n_samples)
    else
        error("Unknown method")
    end

    return maxlambda
end

function _maxlambda4_brute_force(m::AbstractMatrix)
    n = size(m, 2)
    n_include = ceil(Int, n / 2)
    is = axes(m, 2)

    combs = combinations(is, n_include)
    ncombs = length(combs)

    if ncombs > 1e6
        @info "Brute forcing $(ncombs) combinatinos. Adjust your expectations accordingly..."
    end

    maxlambda = maximum(lambda4(m, c) for c in combs)

    return maxlambda
end

function _maxlambda4_random(m::AbstractMatrix, n_samples::Int)
    n = size(m, 2)
    n_include = ceil(Int, n / 2)
    is = axes(m, 2)
    maxlambda =
        maximum(lambda4(m, sample(is, n_include, replace = false)) for _ in 1:n_samples)
    return maxlambda
end

@kwdef struct L4 <: ReliabilityMeasure
    method::Symbol = :auto
    n_samples::Int = 10_000
end

(method::L4)(data) = maxlambda4(data, method = method.method, n_samples = method.n_samples)
name(r::L4) = "L4(:$(r.method), $(r.n_samples))"

"""
    lambda5(m::AbstractMatrix)
    lambda5(test::PsychometricTest)
    lambda5(test::PsychometricTest, scale::Symbol)

Return the lower bound estimate of the reliability lambda₅ described in $GUTTMAN1945
"""
function lambda5(m::AbstractMatrix)
    C = LowerTriangular(cov(m))
    zerodiag!(C)
    sj = sum(abs2, C, dims = 1)
    Cmax = maximum(sj)
    st = var(sum(m, dims = 2))
    return lambda1(m) + 2 * sqrt(Cmax) / st
end

struct L5 <: ReliabilityMeasure end
(method::L5)(data) = lambda5(data)

"""
    lambda6(m::AbstractMatrix)
    lambda6(test::PsychometricTest)
    lambda6(test::PsychometricTest, scale::Symbol)

Return the lower bound estimate of the reliability lambda₆ described in $GUTTMAN1945
"""
function lambda6(m::AbstractMatrix)
    C = cov(m)
    Cinv = inv(C)
    smc = 1 .- 1 ./ diag(Cinv)
    return 1 - sum(1 .- smc) / sum(C)
end

struct L6 <: ReliabilityMeasure end
(method::L6)(data) = lambda6(data)

"""
    GUTTMAN_METHODS

A collection of reliability measures described in $GUTTMAN1945
"""
const GUTTMAN_METHODS = [L1(), L2(), L3(), L4(), L5(), L6()]
