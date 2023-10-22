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
    numthreads = Threads.nthreads()
    println("Threads detected: $numthreads")
    numcores = min(1,numthreads÷2)
    BLAS.set_num_threads(numcores)
    FFTW.set_num_threads(numcores)
end



function copy_odds_branchless!(dst::Vector{T}, src::Vector{T}) where {T <: Integer}
    write_index = 1
    @inbounds for i in eachindex(src)
        v = src[i]
        dst[write_index] = v
        write_index += isodd(v)
    end
    return dst
end;


ifelse()





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

@code_native rfft_plan * Amem_sol

@code_native rfft(Amem_sol)

@code_llvm rfft(Amem_sol)

@code_lowered solve(ogprobjac, Rodas5P())
@code_typed solve(ogprobjac, Rodas5P())

using JET
@report_opt solve(ogprobjac, Rodas5P())

@allocations sol.u[1]

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
dubious_rowtuples = []
ampvals = []
for dir in readdir("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000"; join=true)
    for combodir in readdir(dir; join=true)
        if occursin("Raw", basename(combodir))
            # @info "Processing directory $combodir"
            for DFdir in readdir(combodir; join=true)
                # extract dataframes into a vector
                dfarray = read_csvs_in_directory(DFdir)

                for df in dfarray
                    # testrowtuple = copy.(eachrow((testrawdf[testrawdf.amp .> sum(df[1,Between(:L, :A)]), :])))
                    # append!(dubious_rowtuples, testrowtuple)
                    # append!(ampvals, df.amp)
                    # if any(df.amp .> 100.)
                    #     @info "Found a solution with amplitude greater than 100 in $DFdir"
                    #     append!(ampvals, df.amp)
                    # end
                    for row in eachrow(df)
                        if row.amp > sum(row[Between(:L, :A)])
                            # @info "Found a solution with amplitude greater than the total starting concentrations in $DFdir"
                            @info "Find unrealistic amplitude $(row.amp)"
                            push!(dubious_rowtuples, row)
                        end
                    end
                end
            end
        end
    end
end
dubious_rowtuples
maximum(ampvals)



                # Test amplitude columns and save row if amplitude is greater than the total starting concentrations

testrow = dubious_rowtuples[2]


ogprobjac = make_ODE_problem()


#* remake with new initial conditions and new parameters
new_p = [param for param in testrow[Between(:ka1, :DF)]]
new_u = [ic for ic in testrow[Between(:L, :A)]]
new_u0 = [new_u; zeros(length(ogprobjac.u0) - length(new_u))]

newprob = remake(ogprobjac, p = new_p, u0 = new_u0)

sol = solve(newprob, Rodas5P())

using Plots
plot(sol)

sol = solve_odeprob(newprob, [6, 9, 10, 11, 12, 15, 16])

solve_for_fitness_peramp(newprob, [6, 9, 10, 11, 12, 15, 16])