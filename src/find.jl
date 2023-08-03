"""
    find(m::AbstractMatrix, n::Int; criterion = glb, progress = true)
    find(test::PsychometricTest, n::Int; criterion = glb, progress = true)

Perform an exhaustive search to find the subset of `n` items with maximum reliability.
"""
function find(m::AbstractMatrix, args...; kwargs...)
    is = _find(m, args...; kwargs...)
    return m[:, is]
end

function find(test::PsychometricTest, args...; kwargs...)
    responses = response_matrix(test)
    is = _find(responses, args...; kwargs...)
    return subset(test, :, is)
end

function _find(m::AbstractMatrix, n::Int; criterion = glb, progress = true, kwargs...)
    is = axes(m, 2)
    combs = combinations(is, n)

    optimal_is = zeros(Int, n)
    max_crit = -Inf

    prog = Progress(
        length(combs),
        dt = 0.5,
        barglyphs = BarGlyphs("[=> ]"),
        enabled = progress,
    )

    for (i, c) in enumerate(combs)
        subtest = view(m, :, c)
        crit = criterion(subtest, kwargs...)

        if crit > max_crit
            max_crit = crit
            optimal_is = c
        end

        ProgressMeter.update!(prog, i, showvalues = [(:reliability, max_crit)])
    end

    ProgressMeter.finish!(prog)

    return optimal_is
end
