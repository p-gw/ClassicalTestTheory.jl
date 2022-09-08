abstract type AbstractTest end

function scales(test::AbstractTest) end

struct SingleScaleTest{T<:AbstractScale} <: AbstractTest
    scale::T
end

scales(test::SingleScaleTest) = test.scale

struct MultiScaleTest{T<:AbstractScale} <: AbstractTest
    scales::Vector{T}
end

scales(test::MultiScaleTest) = test.scales

Test(scale::AbstractScale) = SingleScaleTest(scale)
Test(data::AbstractMatrix{Bool}) = SingleScaleTest(DichotomousScale("", data))

function Test(data::AbstractMatrix{Int})
    elements = unique(data)
    if all(elements .∈ Ref(0:1))
        scale = DichotomousScale("", data)
    else
        # scale = OrdinalScale("", data)
    end
    return SingleScaleTest(scale)
end

function Test(scales...)
    allunique(id.(scales)) || throw(ArgumentError("All scale ids must be unique"))
    return Test(collect(scales))
end
