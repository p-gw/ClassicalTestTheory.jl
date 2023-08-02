using ClassicalTestTheory
using PsychometricTests
using DelimitedFiles
using CSV
using DataFrames

data = CSV.read("benchmark/attitude.csv", DataFrame)
m = Matrix(data)
test = PsychometricTest(data)

ia = ClassicalTestTheory.itemanalysis(test)

# λ1
λ1(m)
λ1(test)
λ1(test_scales)
λ1(test_scales, nothing)

@benchmark λ1($m)
@benchmark λ1($test)
@benchmark λ1($test_scales)
@benchmark λ1($test_scales, $:a)
@benchmark λ1($test_scales, $nothing)

@profview [λ1(m) for _ in 1:10_000]
@profview [λ1(test) for _ in 1:10_000]
@profview [λ1(test_scales) for _ in 1:10_000]

@profview_allocs λ1(test_scales, :a) sample_rate = 1

@code_warntype λ1(test)
@code_warntype λ1(test_scales)
@code_warntype λ1(test_scales, :a)
