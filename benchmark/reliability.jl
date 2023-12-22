using ClassicalTestTheory
using PsychometricTests
using DelimitedFiles
using CSV
using DataFrames

data = CSV.read("benchmark/attitude.csv", DataFrame)
m = Matrix(data)
test = PsychometricTest(data)

ia = ClassicalTestTheory.itemanalysis(test)

ClassicalTestTheory.reliability(m, ClassicalTestTheory.PSYCH_METHODS)

# lambda1
lambda1(m)
lambda1(test)
lambda1(test_scales)
lambda1(test_scales, nothing)

@benchmark lambda1($m)
@benchmark lambda1($test)
@benchmark lambda1($test_scales)
@benchmark lambda1($test_scales, $:a)
@benchmark lambda1($test_scales, $nothing)

@profview [lambda1(m) for _ in 1:10_000]
@profview [lambda1(test) for _ in 1:10_000]
@profview [lambda1(test_scales) for _ in 1:10_000]

@profview_allocs lambda1(test_scales, :a) sample_rate = 1

@code_warntype lambda1(test)
@code_warntype lambda1(test_scales)
@code_warntype lambda1(test_scales, :a)

# maxlambda4

# find
find(test, 6)

@benchmark find($test, $6)

@code_warntype find(test, 6)

@profview_allocs find(test, 6, criterion = alpha) sample_rate = 1

# test bootstrapping
using Bootstrap

