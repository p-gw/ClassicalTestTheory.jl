using ClassicalTestTheory
using Distributions
using PsychometricTests
using Test

@testset "ClassicalTestTheory.jl" begin
    include("reliability.jl")
    include("split.jl")
end
