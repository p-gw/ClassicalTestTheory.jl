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
    item_ids = getid.(getitems(test)[is])
    return subset(test, :, item_ids)
end

function find(test::PsychometricTest, scale::Symbol, args...; kwargs...)
    responses = response_matrix(test, scale)
    is = _find(responses, args...; kwargs...)
    item_ids = getid.(getitems(test)[is])
    return subset(test, :, item_ids)
end

function _find(
    m::AbstractMatrix,
    n::Int;
    criterion::F = glb,
    progress = true,
    kwargs...,
) where {F}
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
    max_crit = -Inf

    prog = Progress(
        length(combs),
        dt = 0.5,
        barglyphs = BarGlyphs("[=> ]"),
        enabled = progress,
    )

    for (i, c) in enumerate(combs)
        subtest = view(m, :, c)
        crit = criterion(subtest; kwargs...)

        if crit > max_crit
            max_crit = crit
            optimal_is = c
        end

        ProgressMeter.update!(
            prog,
            i,
            showvalues = [(:items, optimal_is), (:reliability, max_crit)],
        )
    end

    ProgressMeter.finish!(prog)

    return optimal_is
end
