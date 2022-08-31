

"""
    split(m::AbstractMatrix, is)

Split a matrix by arbitrary indices `is`.

Returns a `Tuple` of matrices `m1`, `m2` with `size(m1, dims=1) == length(is)` and `size(m2, dims=1) == size(m, dims=1) - length(is)`.
"""

function Base.split(m::AbstractMatrix, is)
    n = nitems(m)
    exclude = getexcludeindices(n, is)
    return view(m, :, is), view(m, :, exclude)
end

"""
    split(t::Test, is)

Split a test by arbitrary indices `is`.
"""
function Base.split(t::Test, is)
    n = nitems(t)
    exclude = getexcludeindices(n, is)
    t1 = SubTest(t, is)
    t2 = SubTest(t, exclude)
    return t1, t2
end

getexcludeindices(n, is::Union{BitVector,Vector{Bool}}) = .!is
getexcludeindices(n, is) = setdiff(1:n, is)

"""
    splithalf(m::AbstractMatrix; type::Symbol)

Split a matrix rowwise in two halfes.
The type of split is determined by `type`

# Available types
- `:oddeven`: Split the matrix by odd and even indices
- `:firstlast`: Split the matrix by first and last half
- `:random`: Split the matrix by random indices
"""
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

function splithalf(m::AbstractMatrix; kwargs...)
    col_is = 1:size(m, 2)
    is = getsplitindices(col_is; kwargs...)
    return m[:, is], m[:, .!is]
end

"""
"""
function splithalf(t::Test; kwargs...)
    col_is = 1:nitems(t)
    is = getsplitindices(col_is, kwargs...)
    t1 = SubTest(t, is)
    t2 = SubTest(t, .!is)
    return t1, t2
end

_oddeven_is(is) = iseven.(is)
_firstlast_is(is) = is .<= ceil(length(is) / 2)
function _random_is(is)
    s = sample(is, ceil(Int, length(is) / 2), replace=false)
    return [i in s for i in is]
end
