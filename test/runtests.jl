using ClassicalTestTheory
using Distributions
using Test

@testset "ClassicalTestTheory.jl" begin
    include("reliability.jl")
    include("split.jl")
    include("find.jl")
end
