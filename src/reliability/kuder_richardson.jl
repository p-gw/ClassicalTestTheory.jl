"""
    kr20(m::AbstractMatrix)
"""
function kr20(m::AbstractMatrix)
    n = size(m, 2)
    item_facilities = mean(m, dims = 1)
    item_difficulties = 1 .- item_facilities
    st = var(sum(m, dims = 2))
    return (n / (n - 1)) * ((st - sum(item_facilities .* item_difficulties)) / st)
end

struct KR20 <: ReliabilityMeasure end
(method::KR20)(data) = kr20(data)

"""
    kr21(m::AbstractMatrix)
"""
function kr21(m::AbstractMatrix)
    n = size(m, 2)
    item_facilities = mean(m, dims = 1)
    avg_facility = mean(item_facilities)
    avg_difficulty = 1 - avg_facility
    st = var(sum(m, dims = 2))
    return (n / (n - 1)) * ((st - n * avg_difficulty * avg_facility) / st)
end

struct KR21 <: ReliabilityMeasure end
(method::KR21)(data) = kr21(data)
