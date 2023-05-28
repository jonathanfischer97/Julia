using BenchmarkTools, Profile

testarray = rand(Float64, 1000000)
testpeakidxs, testpeakvals = findmaxima(testarray, 1)

@btime @fastmath getDif($testpeakidxs, $testarray)


#< Benchmarking automatic vs stiff solvers
@benchmark autosol = solve(fullprob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
@benchmark @fastmath solve(fullprob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
plot(autosol)
@benchmark stiffsol = solve(fullprob, Rosenbrock23(), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
plot!(stiffsol)

difference = sum(autosol - stiffsol)



@btime CostFunction($autosol)
@btime @fastmath CostFunction($autosol)


#< PerAmp testing 
per, amp = getPerAmp(autosol)


eval_ic_fitness(fullprob.u0, fullprob)


#< Comparing getDif_bidirectional to getDif
getDif(testpeakidxs, testarray)
getDif(testpeakvals)
getDif_bidirectional(testpeakvals)


#< Testing Supertypes 
abstract type AbstractTestType <: ConstraintType end

struct TestType1 <: AbstractTestType 
    a::Int
    b::Int
    c::Int
end

testype1 = TestType1(1,2,3)

for n in testype1
    println(n)
end