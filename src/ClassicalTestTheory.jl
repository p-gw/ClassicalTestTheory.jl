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

export Test, SubTest
export scores, responses, nitems, npersons

export split, splithalf

export difficulty, facility, itc

export λ1, λ2, λ3, λ4, maxλ4, λ5, λ6
export α
export kr20, kr21
export glb

include("utils.jl")
include("types.jl")
include("split.jl")
include("item_statistics.jl")
include("reliability.jl")

end
