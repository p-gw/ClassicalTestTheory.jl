abstract type AbstractTest end

struct Test{T<:Real} <: AbstractTest
    data::Matrix{T}
    itemcov::Matrix{Float64}
    scores::Vector{T}
end

struct SubTest <: AbstractTest
    data
    itemcov
    scores
end

Test(m::AbstractMatrix) = Test{eltype(m)}(m, cov(m), scores(m))
Test(st::SubTest) = Test{eltype(st.data)}(st.data, st.itemcov, st.scores)

function SubTest(t::Test, is)
    dview = view(t.data, :, is)
    return SubTest(
        dview,
        view(t.itemcov, is, is),
        scores(dview)
    )
end

scores(v::AbstractVector) = v
scores(m::AbstractMatrix) = vec(sum(m, dims=2))
scores(t::AbstractTest) = t.scores

responses(m::AbstractMatrix) = m
responses(t::AbstractTest) = t.data
responses(x, i::Int) = responses(x)[:, i]

nitems(m::AbstractMatrix) = size(m, 2)
nitems(t::AbstractTest) = nitems(t.data)

npersons(m::AbstractMatrix) = size(m, 1)
npersons(t::AbstractTest) = npersons(t.data)

Statistics.cov(t::AbstractTest) = t.itemcov
