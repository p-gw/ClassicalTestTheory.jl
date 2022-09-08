const GUTTMAN1945 = "Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255-282."

"""
   λ1(scale::AbstractScale)::Float64
   λ1(test::SingleScaleTest)::Float64
   λ1(test::MultiScaleTest)::Vector{Float64}

Return the lower bound estimate of the reliability L₁ described in $GUTTMAN1945

``\\lambda_1 = 1 - \\frac{\\sum_{j=1}^{n} s_j^2}{s_t^2}``
"""
function λ1(scale::AbstractScale)
    sum_sj = tr(cov(scale))
    st = var(scores(scale))
    return 1 - (sum_sj / st)
end

λ1(test::SingleScaleTest) = λ1(scales(test))
λ1(test::MultiScaleTest) = λ1.(scales(test))

"""
    λ2(scale::AbstractScale)::Float64
    λ2(test::SingleScaleTest)::Float64
    λ2(test::MultiScaleTest)::Vector{Float64}

Return the lower bound estimate of the reliability λ₂ described in $GUTTMAN1945

``\\lambda_2 = \\lambda_1 + \\frac{\\sqrt{\\frac{n}{n-1}C_2}}{s_t^2}``
"""
function λ2(scale::AbstractScale)
    n = nitems(scale)
    C2 = cov(scale) .^ 2
    zerodiag!(C2)
    st = var(scores(scale))
    return λ1(scale) + sqrt(n / (n - 1) * sum(C2)) / st
end

λ2(test::SingleScaleTest) = λ2(scales(test))
λ2(test::MultiScaleTest) = λ2.(scales(test))

"""
    λ3(scale::AbstractScale)::Float64
    λ3(test::SingleScaleTest)::Float64
    λ3(test::MultiScaleTest)::Vector{Float64}

Return the lower bound estimate of the reliability λ₃ described in $GUTTMAN1945

``\\lambda_3 = \\frac{n}{n-1}\\lambda_1``
"""
function λ3(scale::AbstractScale)
    n = nitems(scale)
    return n / (n - 1) * λ1(scale)
end

λ3(test::SingleScaleTest) = λ3(scales(test))
λ3(test::MultiScaleTest) = λ3.(scales(test))

"""
    α(x)

Estimate Cronbach's α. `α` is an alias for [`λ3`](@ref).
"""
const α = λ3

"""
    λ4(scale::AbstractScale; kwargs...)::Float64
    λ4(test::SingleScaleTest; kwargs...)::Float64
    λ4(test::MultiScaleTest; kwargs...)::Vector{Float64}

Return the lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

``\\lambda_4 = 2\\( 1 - \\frac{s_a^2 + s_b^2}{s_t^2} \\)``

The calculation of λ₄ is based on splitting `scale` in half.
It is a lower bound of the reliability no matter how the scale is split.

The split of the scale can be controlled by the `type` keyword argument (see also [`splithalf`](@ref)).

To get the maximum lower bound see [`maxλ4`](@ref).
"""
function λ4(scale::AbstractScale; kwargs...)
    subscales = splithalf(scale; kwargs...)
    st = var(scores(scale))
    vars = var.(scores.(subscales))
    return 2 * (1 - sum(vars) / st)
end

λ4(test::SingleScaleTest; kwargs...) = λ4(scales(test); kwargs)
λ4(test::MultiScaleTest; kwargs...) = λ4.(scales(test); kwargs)

"""
"""
function λ4(scale::AbstractScale, is)
    subtests = split(x, is)
    st = var(scores(x))
    vars = var.(scores.(subtests))
    return 2 * (1 - sum(vars) / st)
end

λ4(test::SingleScaleTest, is) = λ4(scales(test), is)

"""
    maxλ4(x; method=:auto, n_samples=1_000)

Return the maximum lower bound estimate of the reliability λ₄ described in $GUTTMAN1945

By default (if `n_samples=nothing`) the maximum value is found by brute force iteration over
the split-half permutations of `x`.
As the number of permutations grows, this method becomes infeasible.
For large numbers of items a (local) maximum of λ₄ can be found by sampling random split-half
permutations.
To use sampling, specify `n_samples::Int`.

See also [`λ4`](@ref).
"""
function maxλ4(x; method=:auto, n_samples=10_000)::Float64
    if method == :auto
        if nitems(x) <= 25
            maxλ = _maxλ4_brute_force(x)
        else
            maxλ = _maxλ4_random(x, n_samples)
        end
    elseif method == :bruteforce
        maxλ = _maxλ4_brute_force(x)
    elseif method == :sample
        maxλ = _maxλ4_random(x, n_samples)
    else
        error("Unknown method")
    end

    return maxλ
end

function _maxλ4_brute_force(x)
    n = nitems(x)
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
    maxλ = maximum(λ4(x, perm) for perm in Iterators.take(perms, n_perms_iter))
    return maxλ
end

function _maxλ4_random(x, n_samples::Int)
    n = nitems(x)
    n_include = ceil(Int, n / 2)
    is = vcat(trues(n_include), falses(n - n_include))
    maxλ = maximum(λ4(x, shuffle!(is)) for _ in 1:n_samples)
    return maxλ
end

"""
    λ5(x)

Return the lower bound estimate of the reliability λ₅ described in $GUTTMAN1945

``\\lambda_5 = \\lambda_1 + \\frac{2\\sqrt{\\bar{C}_2}}{s_t^2}``
"""
function λ5(x)
    covmat = LowerTriangular(copy(cov(x)))
    zerodiag!(covmat)
    C2 = sum(covmat .^ 2, dims=1)
    C2_max = maximum(C2)
    st = var(scores(x))
    return λ1(x) + 2 * sqrt(C2_max) / st
end

"""
    λ6(x)

Return the lower bound estimate of the reliability λ₆ described in $GUTTMAN1945

``\\lambda_6 = 1 - \\frac{\\sum_{j=1}^n e_j^2}{s_t^2}``
"""
function λ6(x)
    covmat = cov(x)
    inv_covmat = inv(covmat)
    smc = 1 .- 1 ./ diag(inv_covmat)
    return 1 - sum(1 .- smc) / sum(covmat)
end

"""
    kr20
"""
function kr20(x)
    n = nitems(x)
    item_facilities = facility.(eachitem(x))
    item_difficulties = 1 .- item_facilities
    st = var(scores(x))
    return (n / (n - 1)) * ((st - sum(item_facilities .* item_difficulties)) / st)
end

"""
    kr21
"""
function kr21(x)
    n = nitems(x)
    avg_facility = mean(facility.(eachitem(x)))
    avg_difficulty = 1 - avg_facility
    st = var(scores(x))
    return (n / (n - 1)) * ((st - n * avg_difficulty * avg_facility) / st)
end

"""
    glb(x)

Return the greatest lower bound estimate (glb) of the reliability as described in
Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items: II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579-591.
"""
function glb(x)
    n = nitems(x)

    C = cov(x)
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
