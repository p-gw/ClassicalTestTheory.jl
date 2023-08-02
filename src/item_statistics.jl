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
function itc(x, i; corrected = true, standardize = false)
    if corrected in (true, :henrysson)
        return _itc_henrysson(x, i; standardize)
    elseif corrected == false
        return _itc_uncorrected(x, i; standardize)
    elseif corrected == :zubin
        return _itc_zubin(x, i; standardize)
    elseif corrected == :guilford
        return _itc_guilford(x, i; standardize)
    else
        error("Unknown correction type: $(corrected)")
    end
end

function _itc_uncorrected(m::AbstractMatrix, i; standardize)
    responses = vec(m[:, i])
    total_scores = vec(sum(m, dims = 2))
    return cor(responses, total_scores)
end

function _itc_uncorrected(test::PsychometricTest, i; standardize)
    responses = response_matrix(test)
    item_ind = findfirst(x -> x == i, getid.(getitems(test)))
    return _itc_uncorrected(responses, item_ind; standardize)
end

function _itc_henrysson(m::AbstractMatrix, i; standardize)
    is = filter(x -> x != i, axes(m, 2))
    corrected_m = view(m, :, is)
    responses = vec(m[:, i])
    total_scores = vec(sum(corrected_m, dims = 2))
    return cor(responses, total_scores)
end

function _itc_henrysson(test::PsychometricTest, i; standardize)
    responses = response_matrix(test)
    item_ind = findfirst(x -> x == i, getid.(getitems(test)))
    return _itc_henrysson(responses, item_ind; standardize)
end

# function _itc_zubin(x, i)
#     itc = _itc_uncorrected(x, i)
#     st = std(scores(x))
#     fv = facility(x, i)
#     y = cdf(Normal(), fv)
#     pqy = fv * (1 - fv) / y
#     return (itc * st - pqy) / sqrt(st^2 + 1 - 2 * itc * st)
# end

# function _itc_guilford(x, i)
#     itc = _itc_uncorrected(x, i)
#     st = std(scores(x))
#     fv = facility(x, i)
#     y = cdf(Normal(), fv)
#     pqy = fv * (1 - fv) / y
#     return (itc * st - pqy) / sqrt(st^2 + pqy^2 - 2 * itc * st * pqy)
# end

struct ItemStatistics{Ti,T<:Real}
    item::Ti
    n::Int
    missings::Int
    itc::T
    itc_std::T
    itc_corrected::T
    mean::T
    std::T
end

function ItemStatistics(test::PsychometricTest, i)
    item_responses = vec(getvalue.(test[:, i]))

    stats = ItemStatistics(
        i,
        sum(!ismissing, item_responses),  # number of complete cases
        sum(ismissing, item_responses),  # number of missings
        itc(test, i, corrected = false),
        itc(test, i, corrected = false, standardize = true),
        itc(test, i, corrected = :henrysson),
        itemmean(test, i),
        std(item_responses),
    )

    return stats
end

struct ItemAnalysis
    statistics::Vector{ItemStatistics}
end

function Base.show(io::IO, items::ItemAnalysis)
    tbl = Tables.columns(items.statistics)

    # TODO: nicer solution
    tbl.itc .= round.(tbl.itc, digits = 2)
    tbl.itc_std .= round.(tbl.itc_std, digits = 2)
    tbl.itc_corrected .= round.(tbl.itc_corrected, digits = 2)
    tbl.mean .= round.(tbl.mean, digits = 2)
    tbl.std .= round.(tbl.std, digits = 2)

    header = ["item", "N", "missings", "itc", "itc (std)", "itc (cor)", "mean", "std"]
    print(
        io,
        # Term.Panel(
        Term.Table(tbl; header, compact = false),
        # title = "{dim}Item Analysis",
        # style = "dim",
        # subtitle = "{dim}ClassicalTestTheory.jl",
        # subtitle_justify = :right,
        # ),
    )

    return nothing
end

function itemanalysis(test::PsychometricTest)
    item_statistics = [ItemStatistics(test, getid(i)) for i in getitems(test)]
    return ItemAnalysis(item_statistics)
end
