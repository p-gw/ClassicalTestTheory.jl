using ClassicalTestTheory
using Distributions
using Random
using StatsBase
using Statistics
using Test

@testset "ClassicalTestTheory.jl" begin
    include("scale.jl")
    include("test.jl")
    # include("types.jl")
    # include("reliability.jl")
    # include("split.jl")
end
