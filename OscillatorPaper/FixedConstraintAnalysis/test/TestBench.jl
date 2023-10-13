begin 
    # using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    # using DiffEqCallbacks
    using Statistics

    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter

    # using JET

    using LinearAlgebra

    using MultivariateStats
    using Clustering

    # using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)



    using OscTools 


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW
    FFTW.set_num_threads(18)
end

# function test4fixedGA()
#     #* Set up the default GA problem
#     ga_problem = GAProblem()

#     #* Fixed some constraints 
#     set_fixed_constraints!(ga_problem.constraints, [:DF, :K, :P, :A])

#     #* Assign the fixed values 
#     set_fixed_values!(ga_problem.constraints, 1000., 1.0, 1.0, 3.16)

#     #* Set seed 
#     Random.seed!(1234)

#     #* Generate the initial population
#     population = generate_population(ga_problem.constraints, 10000)

#     #* Run the GA
#     run_GA(ga_problem, population)
# end

ga_result = test4fixedGA()

population_to_matrix(ga_result)

@benchmark test4fixedGA(5000)
"""BenchmarkTools.Trial: 1 sample with 1 evaluation.
Single result which took 26.816 s (32.50% GC) to evaluate,
with a memory estimate of 92.33 GiB, over 547954601 allocations."""


#*FFTW testing 

ogprobjac = make_ODE_problem()

sol = solve_odeprob(ogprobjac, [6, 9, 10, 11, 12, 15, 16])
Amem_sol = map(sum, sol.u)

"""Returns RFFT plan for the given ODEProblem, using the default problem solution"""
function make_rfft_plan(ode_problem::OP) where OP <: ODEProblem
    sol = solve_odeprob(ode_problem, [6, 9, 10, 11, 12, 15, 16])
    Amem_sol = map(sum, sol.u)
    return plan_rfft(Amem_sol)
end

plan = make_rfft_plan(ogprobjac)

@btime rfft($Amem_sol)
rfft_plan = plan_rfft(Amem_sol)

@btime rfft_plan * $Amem_sol

@btime dct($Amem_sol)

dct_plan = plan_dct(Amem_sol)

@btime dct_plan * $Amem_sol

dct_plan = plan_dct!(Amem_sol)

@btime dct_plan * $Amem_sol

testrawdf = CSV.read("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000/ka1_K_A/RawData/DF=100.0/ka1=0.0_K=1.0_A=0.56.csv", DataFrame)
testrowtuple = copy.(eachrow((testrawdf[testrawdf.amp .< sum(testrawdf[1,Between(:L, :A)]), :])))
testrowtuple[1]

#< TESTING AMPLITUDE CALCULATIONS ##
#* Extract every raw solution that has an amplitude greater than the total starting concentrations
for dir in readdir("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000"; join=true)
    if occursin("Raw", basename(dir))
        for DFdir in readdir(dir; join=true)
            # extract dataframes into a vector
            dfarray = read_csvs_in_directory(DFdir)

            for df in dfarray
                df[(testrawdf.amp .> sum(testrawdf[1,Between(:L, :A)])), :]


                # Test amplitude columns and save row if amplitude is greater than the total starting concentrations



