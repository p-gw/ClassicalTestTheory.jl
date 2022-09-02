"""
    split(m::AbstractMatrix, is)

Split a matrix by arbitrary indices `is`.

Returns a `Tuple` of matrices `m1`, `m2` with `size(m1, dims=1) == length(is)` and `size(m2, dims=1) == size(m, dims=1) - length(is)`.
"""
function Base.split(m::AbstractMatrix, is)
    m1 = view(m, :, is)
    m2 = view(m, :, Not(is))
    return m1, m2
end

"""
    split(t::Test, is)

Split a test by arbitrary indices `is`.
"""
function Base.split(t::AbstractTest, is)
    t1 = SubTest(t, is)
    t2 = SubTest(t, Not(is))
    return t1, t2
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
function splithalf(x; kwargs...)
    n = nitems(x)
    is = getsplitindices(1:n; kwargs...)
    return split(x, is)
end

function getsplitindices(x; type::Symbol=:firstlast)
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

_oddeven_is(is) = iseven.(is)
_firstlast_is(is) = is .<= ceil(length(is) / 2)

function _random_is(is)
    s = sample(is, ceil(Int, length(is) / 2), replace=false)
    return [i in s for i in is]
end
