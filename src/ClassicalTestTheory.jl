module ClassicalTestTheory

using Base: split
using Combinatorics
using SCS
using Distributions
using InvertedIndices
using JuMP
using LinearAlgebra
using Random
using StatsBase
using Statistics
using PsychometricTests
using ProgressMeter
using Memoization
using Term
using Tables
using Printf
using ThreadsX
using SplittablesBase

export Test, SubTest
export eachitem, eachperson
export scores, responses, nitems, npersons

export split, splithalf

export difficulty, facility, itc

export λ1, mλ1, λ2, λ3, λ4, maxλ4, λ5, λ6
export α
export kr20, kr21
export glb
export find

# new
export DichotomousScale, OrdinalScale
export id, responses, scores

export Test
export scales

include("references.jl")
include("utils.jl")
include("reliability.jl")

# include("types.jl")

# include("scale.jl")
# include("test.jl")

include("split.jl")

include("item_statistics.jl")
include("find.jl")

end
