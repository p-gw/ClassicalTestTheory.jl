"""
    find(m::AbstractMatrix, n::Int, method::ReliabilityMeasure = GLB(); progress = true)

Perform an exhaustive search to find the subset of `n` items with maximum reliability, where
`method` is used to estimate the reliability.
"""
function find(m::AbstractMatrix, args...; kwargs...)
    is = _find(m, args...; kwargs...)
    return m[:, is]
end

function _find(
    m::AbstractMatrix,
    n::Int,
    method::ReliabilityMeasure = GLB();
    progress = true,
)
    if n >= size(m, 2)
        throw(
            ArgumentError(
                "The subset size must be smaller than the size of the orginial test.",
            ),
        )
    end

    is = axes(m, 2)
    combs = combinations(is, n)

    optimal_is = zeros(Int, n)
    max_reliability = Ref(-Inf)

    if progress
        prog = ProgressBar(transient = true)

        Progress.with(prog) do
            prog_job = addjob!(
                prog,
                N = length(combs),
                description = "Finding optimal item subset...",
            )

            for c in combs
                _update_reliability!(max_reliability, optimal_is, m, c, method)
                update!(prog_job)
            end
        end
    else
        for c in combs
            _update_reliability!(max_reliability, optimal_is, m, c, method)
        end
    end

    return optimal_is
end

function _update_reliability!(max_reliability, optimal_is, m, c, method)
    subtest = view(m, :, c)
    reliability = method(subtest)

    if reliability > max_reliability[]
        max_reliability[] = reliability
        optimal_is .= c
    end

    return nothing
end
