module ClassicalTestTheory

using Base: @kwdef
using Bootstrap
using Combinatorics
using Distributions
using JuMP
using LinearAlgebra
using OrderedCollections
using Printf
using ProgressMeter
using PsychometricTests
using Random
using Reexport
using SCS
using Statistics
using StatsAPI
using StatsBase
using Tables
using Term

@reexport import StatsAPI: confint, stderror
import Base: split

# reliability measures
export lambda1, lambda2, lambda3, lambda4, maxlambda4, lambda5, lambda6, alpha
export L1, L2, L3, L4, L5, L6, Alpha
export GUTTMAN_METHODS, PSYCH_METHODS, mu_up_to

export kr20, kr21
export KR20, KR21

export glb, mu
export GLB, Mu

export reliability

export estimate, bootstrap_sample

# item statistics
export itc, itemanalysis

# find
export find

include("references.jl")
include("utils.jl")
include("split.jl")
include("reliability/reliability.jl")
include("item_statistics.jl")
include("find.jl")

end
