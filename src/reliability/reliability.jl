const DEFAULT_BOOSTRAP_CI = BCaConfInt(0.95)
const DEFAULT_BOOSTRAP_SAMPLING = BasicSampling(1000)

"""
    ReliabilityMeasure

A reliability measure is the reliability interface for end users.
All low-level reliability estimates should define a callable struct `T <: ReliabilityMeasure`
that takes a single argument, the data, and return the reliability estimate.

## Methods
- `(method::T)(data)`: The reliability algorithm
- `name(method::T)`: The name of the reliability algorithm for pretty printing
"""
abstract type ReliabilityMeasure <: Function end

name(r::ReliabilityMeasure) = string(nameof(r))

include("guttman.jl")
include("glb.jl")
include("kuder_richardson.jl")
include("ten_berge.jl")

"""
    PSYCH_METHODS

A collection of reliability measures used by the `psych` R package.
"""
const PSYCH_METHODS = [L4(), L6(), L3(), L2()]

"""
    ReliabilityResult

A data structure that holds results from reliability analysis.

## Methods
- `bias`
- `ci_method`
- `confint`
- `data`
- `estimate`
- `samples`
- `sampling_method`
- `stderror`
"""
struct ReliabilityResult{D,B,S,C}
    data::D
    bootstrap_samples::B
    sampling_method::S
    ci_method::C
end

"""
    data(result::ReliabilityResult)

Get the data used to calculate `result`.
"""
data(result::ReliabilityResult) = result.data

"""
    sampling_method(result::ReliabilityResult)

Get the sampling method used to calculate the bootstrap samples in `result`.
"""
sampling_method(result::ReliabilityResult) = result.sampling_method

"""
    ci_method(result::ReliabilityResult)

Get the method used to calculate the bootstrap confidence interval in `result`.
"""
ci_method(result::ReliabilityResult) = result.ci_method

"""
    estimate(result::ReliabilityResult)::OrderedDict{String, Float64}
    estimate(result::ReliabilityResult, method::String)::Float64

Get the point estimate of the reliability.
"""
function estimate(result::ReliabilityResult, method::String)
    samples = result.bootstrap_samples[method]
    return first(original(samples))
end

function estimate(result::ReliabilityResult)
    samples = result.bootstrap_samples
    return OrderedDict(key => estimate(result, key) for key in keys(samples))
end

"""
    confint(result::ReliabilityResult)::OrderedDict{String, Tuple{Float64, Float64}}
    confint(result::ReliabilityResult, method::String)::Float64

Get the bootstrap confidence intervals for the reliability estimates of `result`.
"""
function StatsAPI.confint(result::ReliabilityResult, method::String)
    samples = result.bootstrap_samples
    ci = confint(samples[method], result.ci_method)
    _, lwr, upr = first(ci)
    return lwr, upr
end

function StatsAPI.confint(result::ReliabilityResult)
    samples = result.bootstrap_samples
    return OrderedDict(key => confint(result, key) for key in keys(samples))
end

"""
    bias
"""
function Bootstrap.bias(r::ReliabilityResult, method::String)
    b = Bootstrap.bias(r.bootstrap_samples[method])
    return first(b)
end

function Bootstrap.bias(r::ReliabilityResult)
    return OrderedDict(key => bias(r, key) for key in keys(r.bootstrap_samples))
end

"""
    stderror
"""
function Bootstrap.stderror(r::ReliabilityResult, method::String)
    b = Bootstrap.stderror(r.bootstrap_samples[method])
    return first(b)
end

function Bootstrap.stderror(r::ReliabilityResult)
    return OrderedDict(key => stderror(r, key) for key in keys(r.bootstrap_samples))
end

"""
    samples
"""
function samples(r::ReliabilityResult, method::String)
    b = straps(r.bootstrap_samples[method])
    return first(b)
end

function samples(r::ReliabilityResult)
    return OrderedDict(key => samples(r, key) for key in keys(r.bootstrap_samples))
end

"""
    reliability(m::AbstractMatrix, method::ReliabilityMeasure; kwargs...)
    reliability(m::AbstractMatrix, methods::Vector{<:ReliabilityMeasure}; kwargs...)

Estimate the reliability of `m` for a given `method` or multiple `methods`.

## Arguments
- `m`: The input data
- `method`: A `ReliabilityMeasure` to estimate

## Keyword arguments
- `ci_method`: The method used to calculate the bootstrap confidence intervals. Defaults to `$(DEFAULT_BOOSTRAP_CI)`
- `sampling_method`: The method used to draw boostrap samples from `m`. Defaults to `$(DEFAULT_BOOSTRAP_SAMPLING)`
"""
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

function reliability(m::AbstractMatrix, method::ReliabilityMeasure; kwargs...)
    return reliability(m, [method]; kwargs...)
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
    return nothing
end
