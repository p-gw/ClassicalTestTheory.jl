const DEFAULT_BOOSTRAP_CI = BCaConfInt(0.95)
const DEFAULT_BOOSTRAP_SAMPLING = BasicSampling(1000)

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

"""
    alpha(m::AbstractMatrix)
    alpha(test::PsychometricTest)
    alpha(test::PsychometricTest, scale::Symbol)

Estimate Cronbach's alpha. `alpha` is an alias for [`lambda3`](@ref).
"""
const alpha = lambda3

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

    model = JuMP.Model(SCS.Optimizer)

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
        error("something went wrong")  # TODO: Fix error message
    end
end

"""
    mu(m::AbstractMatrix, r::Int)
    mu(test::PsychometricTest, r::Int)
    mu(test::PsychometricTest, scale::Symbol, r::Int)

Calculate the lower bound of the reliability mu derived in $TENBERGE1978

## Notes
- If `r = 0` then mu is equivalent to Cronbach's alpha.
- If `r = 1` then mu is equivalent to Guttman's lambda₂.
"""
function mu(m::AbstractMatrix, r::Int)
    r >= 0 || throw(ArgumentError("r must be non-negative."))

    n = size(m, 2)
    C = cov(m)
    zerodiag!(C)

    st = var(sum(m, dims = 2))

    p_sum = zero(st)

    for h in Iterators.reverse(0:r)
        p_h = sum(c -> c^(2.0^h), C)

        if h == r
            p_h *= n / (n - 1)
        end

        if h == 0
            p_sum = p_sum + p_h
        else
            p_sum = sqrt(p_sum + p_h)
        end
    end

    return p_sum / st
end

# Reliability Interface for end users
# 1) single function call to reliability(...)
# 2) Should work with a single method, e.g. lambda1
# 3) Should work with multiple methods,
# 4) Sensible defaults while still allowing customization
abstract type ReliabilityMeasure <: Function end

name(r::ReliabilityMeasure) = string(nameof(r))

struct L1 <: ReliabilityMeasure end
(method::L1)(data) = lambda1(data)

struct L2 <: ReliabilityMeasure end
(method::L2)(data) = lambda2(data)

struct L3 <: ReliabilityMeasure end
(method::L3)(data) = lambda3(data)

@kwdef struct L4 <: ReliabilityMeasure
    method::Symbol = :auto
    n_samples::Int = 10_000
end

(method::L4)(data) = maxlambda4(data, method = method.method, n_samples = method.n_samples)
name(r::L4) = "L4($(r.method), $(r.n_samples))"

struct L5 <: ReliabilityMeasure end
(method::L5)(data) = lambda5(data)

struct L6 <: ReliabilityMeasure end
(method::L6)(data) = lambda6(data)

struct KR20 <: ReliabilityMeasure end
(method::KR20)(data) = kr20(data)

struct KR21 <: ReliabilityMeasure end
(method::KR21)(data) = kr21(data)

struct GLB <: ReliabilityMeasure end
(method::GLB)(data) = glb(data)

@kwdef struct Mu <: ReliabilityMeasure
    r::Int
end

(method::Mu)(data) = mu(data, method.r)
name(r::Mu) = "Mu($(r.r))"

# convenience
const GUTTMAN_METHODS = [L1(), L2(), L3(), L4(), L5(), L6()]
const PSYCH_METHODS = [L4(), L6(), L3(), L2()]
mu_up_to(r::Int = 10) = [Mu(i) for i in 0:r]

struct ReliabilityResult{D,B,S<:Bootstrap.BootstrapSampling,C<:Bootstrap.ConfIntMethod}
    data::D
    bootstrap_samples::B
    sampling_method::S
    ci_method::C
end


function reliability(
    m::AbstractMatrix,
    methods::Vector{<:ReliabilityMeasure};
    ci_method = DEFAULT_BOOSTRAP_CI,
    sampling_method = DEFAULT_BOOSTRAP_SAMPLING,
)
    bootstrap_samples = OrderedDict(
        name(method) => bootstrap(method, m, sampling_method) for method in methods
    )

    return ReliabilityResult(m, bootstrap_samples, sampling_method, ci_method)
end

reliability(m, method; kwargs...) = reliability(m, [method]; kwargs...)

# accessors
data(r::ReliabilityResult) = r.data
sampling_method(r::ReliabilityResult) = r.sampling_method
ci_method(r::ReliabilityResult) = r.ci_method

function estimate(r::ReliabilityResult, method::String)
    return first(original(r.bootstrap_samples[method]))
end

function estimate(r::ReliabilityResult)
    return OrderedDict(key => estimate(r, key) for key in keys(r.bootstrap_samples))
end

function StatsAPI.confint(r::ReliabilityResult, method::String)
    ci = confint(r.bootstrap_samples[method], r.ci_method)
    _, lwr, upr = first(ci)
    return lwr, upr
end

function StatsAPI.confint(r::ReliabilityResult)
    return OrderedDict(key => confint(r, key) for key in keys(r.bootstrap_samples))
end

function Bootstrap.bias(r::ReliabilityResult, method::String)
    b = Bootstrap.bias(r.bootstrap_samples[method])
    return first(b)
end

function Bootstrap.bias(r::ReliabilityResult)
    return OrderedDict(key => bias(r, key) for key in keys(r.bootstrap_samples))
end

function Bootstrap.stderror(r::ReliabilityResult, method::String)
    b = Bootstrap.stderror(r.bootstrap_samples[method])
    return first(b)
end

function Bootstrap.stderror(r::ReliabilityResult)
    return OrderedDict(key => stderror(r, key) for key in keys(r.bootstrap_samples))
end

function samples(r::ReliabilityResult, method::String)
    b = straps(r.bootstrap_samples[method])
    return first(b)
end

function samples(r::ReliabilityResult)
    return OrderedDict(key => samples(r, key) for key in keys(r.bootstrap_samples))
end


function Base.show(io::IO, mime::MIME"text/plain", result::ReliabilityResult)
    ci = confint(result)

    fmt = x -> @sprintf "%.2f" x

    ci_method_name = methodname(ci_method(result))
    ci_level = level(ci_method(result))
    quant_lwr, quant_upr = quantilesfromlevel(ci_level)

    sampling_method_name = methodname(sampling_method(result))
    n_samples = nrun(sampling_method(result))

    partable = OrderedDict(
        "method" => collect(keys(result.bootstrap_samples)),
        "estimate" => fmt.(values(estimate(result))),
        "stderror" => fmt.(values(stderror(result))),
        prettyquantile(quant_lwr) => fmt.(first.(values(ci))),
        prettyquantile(quant_upr) => fmt.(last.(values(ci))),
    )

    header = OrderedDict(
        "1" => [
            "confidence interval method:",
            "confidence level:",
            "bootstrap sampling method:",
            "bootstrap iterations",
        ],
        "2" => [
            "{cyan}$(ci_method_name){/cyan}",
            "{magenta}$(ci_level){/magenta}",
            "{cyan}$(sampling_method_name){/cyan}",
            "{magenta}$(n_samples){/magenta}",
        ],
    )

    print(
        io,
        CTTPanel(
            "",
            Term.Table(
                header,
                box = :NONE,
                compact = true,
                show_header = false,
                columns_justify = :left,
            ),
            # "confidence interval method: {cyan}$(ci_method_name){/cyan}",
            # "confidence level: {magenta}$(ci_level){/magenta}",
            # "bootstrap sampling method: {cyan}$(sampling_method_name){/cyan}",
            # "bootstrap iterations: {magenta}$(n_samples){/magenta}",
            Term.Table(
                partable,
                header_style = "green bold",
                columns_style = ["bold red", "", "", "", ""],
                style = "dim",
                box = :SIMPLE,
            ),
            title = "Reliability Analysis",
        ),
    )
    return ""
end

# generate function definitions for PsychometricTest
fs = [
    :lambda1,
    :lambda2,
    :lambda3,
    :lambda4,
    :maxlambda4,
    :lambda5,
    :lambda6,
    :kr20,
    :kr21,
    :glb,
    :mu,
]

for f in fs
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
