"""
    find(m::AbstractMatrix, n::Int, method::ReliabilityMeasure = GLB(); progress = true)

Perform an exhaustive search to find the subset of `n` items with maximum reliability.
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
    max_reliability = -Inf

    prog = Progress(
        length(combs),
        dt = 0.5,
        barglyphs = BarGlyphs("[=> ]"),
        enabled = progress,
    )

    for (i, c) in enumerate(combs)
        subtest = view(m, :, c)
        reliability = method(subtest)

        if reliability > max_reliability
            max_reliability = reliability
            optimal_is = c
        end

        ProgressMeter.update!(
            prog,
            i,
            showvalues = [(:items, optimal_is), (:reliability, max_reliability)],
        )
    end

    ProgressMeter.finish!(prog)

    return optimal_is
end
