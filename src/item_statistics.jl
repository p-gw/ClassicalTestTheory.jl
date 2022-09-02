"""
    facility(x, i::Int)

Return the facility value (proportion correct) of item `i`.

"""
facility(x, i::Int) = mean(responses(x, i))

"""
    facility(x::AbstractVector)

Return the facility of a vector of item responses `x`.
"""
facility(x::AbstractVector) = mean(x)

"""
    difficulty(x, i::Int)

Return the difficulty (proportion incorrect) of item `i`.
"""
difficulty(x, i::Int) = 1 - facility(x, i)

"""
    difficulty(x::AbstractVector)

Return the difficulty of a vector of item responses `x`.
"""
difficulty(x) = 1 - facility(x)

"""
    itc(x, i::Int; correction=true)

Return the item total correlation of item `i`.

# Correction types
- `false`: Calculate the correlation of item `i` and `x`
- `true`: Apply `:henrysson` correction
- `:henrysson`: Calculate the correlation of item `i` and `x` without `i`
- `:zubin`: Apply the correction of Zubin (1934)
- `:guilford`: Apply the correction of Guilford (1954)

# References
Guilford, J. P. (1954). Psychometric methods. New York: McGraw-Hill.

Henrysson, S. (1963). Correction of item-total correlations in item analysis. *Psychometrika, 28*(2), 211-218.

Zubin, J. The method of internal consistency for selecting test items. *J. educ. Psychol,
1934, 25*, 345-356.
"""
function itc(x, i::Int; correction=true)
    if correction in (true, :henrysson)
        return _itc_henrysson(x, i)
    elseif correction == false
        return _itc_uncorrected(x, i)
    elseif correction == :zubin
        return _itc_zubin(x, i)
    elseif correction == :guilford
        return _itc_guilford(x, i)
    else
        error("Unknown correction type: $(correction)")
    end
end

function _itc_uncorrected(x, i)
    return cor(responses(x, i), scores(x))
end

function _itc_henrysson(x, i)
    # here the second value returns the test with all items except i
    _, corrected_test = split(x, i)
    return cor(responses(x, i), scores(corrected_test))
end

function _itc_zubin(x, i)
    itc = _itc_uncorrected(x, i)
    st = std(scores(x))
    fv = facility(x, i)
    y = cdf(Normal(), fv)
    pqy = fv * (1 - fv) / y
    return (itc * st - pqy) / sqrt(st^2 + 1 - 2 * itc * st)
end

function _itc_guilford(x, i)
    itc = _itc_uncorrected(x, i)
    st = std(scores(x))
    fv = facility(x, i)
    y = cdf(Normal(), fv)
    pqy = fv * (1 - fv) / y
    return (itc * st - pqy) / sqrt(st^2 + pqy^2 - 2 * itc * st * pqy)
end
