# ClassicalTestTheory.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://p-gw.github.io/ClassicalTestTheory.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://p-gw.github.io/ClassicalTestTheory.jl/dev)
[![Build Status](https://github.com/p-gw/ClassicalTestTheory.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/p-gw/ClassicalTestTheory.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/p-gw/ClassicalTestTheory.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/p-gw/ClassicalTestTheory.jl)

ClassicalTestTheory.jl is a Julia package for data analysis using [Classical Test Theory](https://en.wikipedia.org/wiki/Classical_test_theory#:~:text=It%20is%20a%20theory%20of,the%20reliability%20of%20psychological%20tests.).

## Installation
```julia
] add https://github.com/p-gw/ClassicalTestTheory.jl.git
```

## Getting started
ClassicalTestTheory.jl provides two entry points to doing data analsis. 
The input data can either be a numeric `Matrix` or a [`PsychometricTest`](https://github.com/JuliaPsychometrics/PsychometricTests.jl). 
While `Matrix` methods provide full functionality, `PsychometricTest` methods provide some 
additional convenience such as scale analysis. 
For details on how to use ClassicalTestTheory.jl with `PsychometricTest` see [XXX](#).

Consider some input data `x`,

```julia-repl
julia> n_persons = 100;
julia> n_items = 8;

julia> x = rand(0:100, n_persons, n_items);
```

we can get some descriptive analysis of the items,

```julia
itemanalysis(x)
```

or estimate the internal consistency (e.g. using Cronbach's alpha)

```julia
reliability(x, Alpha())
```

The package will automatically calculate the coefficient from the data and construct appropriate confidence intervals.

To get multiple estimates of reliability just pass a vector of methods:

```julia
coefficients = [Alpha(), GLB(), Mu(2)]
reliability(x, coefficients)
```
