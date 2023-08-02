const GUTTMAN1945 = "Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255-282."

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

function λ1(test::PsychometricTest)
    scales = getscales(test)

    if length(scales) > 0
        rel = [s => λ1(test, s) for s in keys(scales)]
    else
        rel = λ1(test, nothing)
    end

    return rel
end

function λ1(test::PsychometricTest, scale::Symbol)
    responses = response_matrix(test, scale)
    rel = λ1(responses)
    return rel
end

function λ1(test::PsychometricTest, ::Nothing)
    responses = response_matrix(test)
    rel = λ1(responses)
    return rel
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

function λ2(test::PsychometricTest)
    scales = getscales(test)

    if length(scales) > 0
        rel = [s => λ2(test, s) for s in keys(scales)]
    else
        rel = λ2(test, nothing)
    end

    return rel
end

function λ2(test::PsychometricTest, scale::Symbol)
    responses = response_matrix(test, scale)
    rel = λ2(responses)
    return rel
end

function λ2(test::PsychometricTest, ::Nothing)
    responses = response_matrix(test)
    rel = λ2(responses)
    return rel
end

"""
    λ3(m::AbstractMatrix)
    λ3(test::PsychometricTest)
    λ3(test::PsychometricTest, scale::Symbol)

Return the lower bound estimate of the reliability λ₃ described in $GUTTMAN1945
"""
function λ3(m::AbstractMatrix)
    n = size(m, 2)
    return n / (n - 1) * λ1(m)
end

function λ3(test::PsychometricTest)
    scales = getscales(test)

    if length(scales) > 0
        rel = [s => λ3(test, s) for s in keys(scales)]
    else
        rel = λ3(test, nothing)
    end

    return rel
end

function λ3(test::PsychometricTest, scale::Symbol)
    responses = response_matrix(test, scale)
    rel = λ3(responses)
    return rel
end

function λ3(test::PsychometricTest, ::Nothing)
    responses = response_matrix(test)
    rel = λ3(responses)
    return rel
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

Return the lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

The calculation of λ₄ is based on splitting `scale` in half.
It is a lower bound of the reliability no matter how the scale is split.

The split of the scale can be controlled by the `type` keyword argument.

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

function λ4(test::PsychometricTest; kwargs...)
    scales = getscales(test)

    if length(scales) > 0
        rel = [s => λ4(test, s; kwargs...) for s in keys(scales)]
    else
        rel = λ4(test, nothing)
    end

    return rel
end

function λ4(test::PsychometricTest, scale::Symbol; kwargs...)
    responses = response_matrix(test, scale)
    rel = λ4(responses; kwargs...)
    return rel
end

function λ4(test::PsychometricTest, ::Nothing; kwargs...)
    responses = response_matrix(test)
    rel = λ4(responses; kwargs...)
    return rel
end

"""
    maxλ4(x; method=:auto, n_samples=1_000)

Return the maximum lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

By default (if `n_samples=nothing`) the maximum value is found by brute force iteration over
the split-half permutations of `x`.
As the number of permutations grows, this method becomes infeasible.
For large numbers of items a (local) maximum of λ₄ can be found by sampling random split-half
permutations.
To use sampling, specify `n_samples`.

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
    is = vcat(trues(n_include), falses(n - n_include))

    perms = multiset_permutations(is, length(is))
    n_perms_iter = Int(length(perms) / 2)

    if n_perms_iter > 1e6
        @info "Brute forcing $(n_perms_iter) permutations. Adjust your expectations accordingly..."
    end

    # we only need to iterate over the first half of the permutations, because
    # multiset_permutations(...) returns a sorted vector and λ4 is symmetric with regards
    # to the splits, e.g. [0, 1] and [1, 0] yield identical values of λ4.
    maxλ = maximum(λ4(m, findall(perm)) for perm in Iterators.take(perms, n_perms_iter))
    return maxλ
end

function _maxλ4_random(m::AbstractMatrix, n_samples::Int)
    n = size(m, 2)
    n_include = ceil(Int, n / 2)
    is = vcat(trues(n_include), falses(n - n_include))
    maxλ = maximum(λ4(m, findall(shuffle!(is))) for _ in 1:n_samples)
    return maxλ
end

"""
    λ5(m::AbstractMatrix)

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
    λ6(x)

Return the lower bound estimate of the reliability λ₆ described in $GUTTMAN1945

``\\lambda_6 = 1 - \\frac{\\sum_{j=1}^n e_j^2}{s_t^2}``
"""
function λ6(m::AbstractMatrix)
    C = cov(m)
    Cinv = inv(C)
    smc = 1 .- 1 ./ diag(Cinv)
    return 1 - sum(1 .- smc) / sum(C)
end

"""
    kr20
"""
function kr20(m::AbstractMatrix)
    n = size(m, 2)
    item_facilities = mean(m, dims = 1)
    item_difficulties = 1 .- item_facilities
    st = var(sum(m, dims = 2))
    return (n / (n - 1)) * ((st - sum(item_facilities .* item_difficulties)) / st)
end

"""
    kr21
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
    glb(x)

Return the greatest lower bound estimate (glb) of the reliability as described in
Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items: II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579-591.
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
