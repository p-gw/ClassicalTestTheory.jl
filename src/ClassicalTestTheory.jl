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
using ParallelUtilities
using MultiFloats

export Test, SubTest
export eachitem, eachperson
export scores, responses, nitems, npersons

export split, splithalf

export difficulty, facility, itc


# new
export DichotomousScale, OrdinalScale
export id, responses, scores

export Test
export scales

# reliability
export lambda1, mlambda1, lambda2, lambda3, lambda4, maxlambda4, lambda5, lambda6
export alpha, mu
export kr20, kr21
export glb
export find

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
