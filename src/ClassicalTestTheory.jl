module ClassicalTestTheory

using Base: split
using Combinatorics
using SCS
using Distributions
using JuMP
using LinearAlgebra
using Random
using StatsBase
using Statistics
using PsychometricTests
using ProgressMeter
using Term
using Tables

# reliability measures
export lambda1, lambda2, lambda3, lambda4, maxlambda4, lambda5, lambda6
export kr20, kr21
export alpha, glb, mu

# item statistics
export itc, itemanalysis

# find
export find

include("references.jl")
include("utils.jl")
include("split.jl")
include("reliability.jl")
include("item_statistics.jl")
include("find.jl")

end
