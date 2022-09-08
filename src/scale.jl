abstract type AbstractScale end

id(scale::AbstractScale) = scale.id
responses(scale::AbstractScale) = scale.responses
scores(scale::AbstractScale) = scale.scores
Statistics.cov(scale::AbstractScale) = scale.itemcov

nitems(scale::AbstractScale) = size(responses(scale), 2)


struct DichotomousScale{T<:AbstractArray{Bool}} <: AbstractScale
    id::String
    responses::T
    itemcov::Matrix{Float64}
    scores::Vector{Int}
end

function DichotomousScale(id, responses::AbstractMatrix)
    itemcov = cov(responses)
    scores = vec(sum(responses, dims=2))
    converted_responses = BitArray(responses)
    return DichotomousScale{typeof(converted_responses)}(id, converted_responses, itemcov, scores)
end

function DichotomousScale(id, data, is)
    responses = getindex(data, :, is)

    # convert vector to matrix for single indices
    if length(is) == 1
        responses = reshape(responses, :, 1)
    end

    return DichotomousScale(id, Matrix{Bool}(responses))
end

struct OrdinalScale{T} <: AbstractScale
    id::String
    responses::T
    itemcov::Matrix{Float64}
    scores::Vector{Int}
end
