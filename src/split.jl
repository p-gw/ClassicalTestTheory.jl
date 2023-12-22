"""
    split(m::AbstractMatrix, is)

Split a matrix by arbitrary indices `is`.

Returns a `Tuple` of matrices `m1`, `m2` with `size(m1, dims=1) == length(is)` and `size(m2, dims=1) == size(m, dims=1) - length(is)`.
"""
function Base.split(m::AbstractMatrix, is)
    not_is = setdiff(axes(m, 2), is)
    m1 = view(m, :, is)
    m2 = view(m, :, not_is)
    return m1, m2
end

"""
    splithalf(m::AbstractMatrix; type::Symbol)

Split a matrix rowwise in two halfes.
The type of split is determined by `type`

# Available types
- `:oddeven`: Split the matrix by odd and even indices
- `:firstlast`: Split the matrix by first and last half
- `:random`: Split the matrix by random indices
"""
function splithalf(m::AbstractMatrix; kwargs...)
    n = size(m, 2)
    is = getsplitindices(1:n; kwargs...)
    return split(m, is)
end

function getsplitindices(x; type::Symbol)
    if type == :oddeven
        is = _oddeven_is(x)
    elseif type == :firstlast
        is = _firstlast_is(x)
    elseif type == :random
        is = _random_is(x)
    else
        error("Unknown sampling type: '$type'")
    end
    return is
end

_oddeven_is(is) = filter(iseven, is)
_firstlast_is(is) = filter(i -> i <= ceil(length(is) / 2), is)

function _random_is(is)
    s = sample(is, ceil(Int, length(is) / 2), replace = false, ordered = true)
    return is[s]
end
