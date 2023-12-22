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

function _itc_henrysson(m::AbstractMatrix, i; standardize)
    is = filter(x -> x != i, axes(m, 2))
    corrected_m = view(m, :, is)
    responses = vec(m[:, i])
    total_scores = vec(sum(corrected_m, dims = 2))
    return cor(responses, total_scores)
end

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

function ItemStatistics(name, test_data, item_data, item_index)
    return ItemStatistics(
        string(name),
        sum(!ismissing, item_data),
        sum(ismissing, item_data),
        itc(test_data, item_index, corrected = false),
        itc(test_data, item_index, corrected = false, standardize = true),
        itc(test_data, item_index, corrected = :henrysson),
        mean(item_data),
        std(item_data),
    )
end

function ItemStatistics(m::AbstractMatrix, j; name = "Item " * string(j))
    item_responses = vec(m[:, j])
    return ItemStatistics(name, m, item_responses, j)
end

struct ItemAnalysis
    statistics::Vector{ItemStatistics}
end

function Base.show(io::IO, items::ItemAnalysis)
    tbl = Tables.columns(items.statistics)

    rounded_cols = [:itc, :itc_std, :itc_corrected, :mean, :std]

    for col in rounded_cols
        coldata = Tables.getcolumn(tbl, col)
        coldata .= round.(coldata, digits = 2)  # TODO: keep trailing digits
    end

    header = ["item", "N", "missings", "itc", "itc (std)", "itc (cor)", "mean", "std"]
    print(
        io,
        Term.Panel(
            Term.Table(
                tbl;
                header,
                header_style = "green bold",
                columns_style = ["bold yellow", "", "", "", "", "", "", ""],
                style = "dim",
                box = :SIMPLE,
                compact = false,
            ),
            title = "{dim}Item Analysis",
            style = "dim",
            subtitle = "{dim}ClassicalTestTheory.jl",
            subtitle_justify = :right,
            fit = true,
        ),
    )

    return nothing
end

function itemanalysis(m::AbstractMatrix)
    item_statistics = [ItemStatistics(m, j) for j in axes(m, 2)]
    return ItemAnalysis(item_statistics)
end
