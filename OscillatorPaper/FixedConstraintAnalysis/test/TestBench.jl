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





ga_df = test4fixedGA()

population_to_matrix(ga_result)

#< amp testing
ampvals = ga_df.amp
minamp = minimum(ampvals)
maxamp = maximum(ampvals)
argmin(ga_df.amp)
argmax(ga_df.amp)

plotboth(ga_df[31,:], ogprobjac)

#> end

#< diversity to population size testing
allconstraints = AllConstraints()
test_population = generate_population(allconstraints, 10000)
testpop = stack(test_population)

# Apply PCA to reduce dimensions
function applyPCA(dataMatrix)
    # pca_model = fit(PCA, dfmat, maxoutdim=3)
    # reduced_data = MultivariateStats.transform(pca_model, dfmat)
    # Center the data
    logdataMatrix = log.(dataMatrix)
    # logdataMatrix = dataMatrix
    centeredData = dataMatrix .- mean(logdataMatrix, dims=1)
    
    # Calculate covariance matrix and its eigen decomposition
    covMatrix = cov(centeredData, dims=1)
    eigenDecomp = eigen(covMatrix)
    
    # Project data into lower-dimensional space
    return centeredData * eigenDecomp.vectors
end

# Calculate the average range of the reduced dimensions
function calculateSpread(reducedData)
    mins = minimum(reducedData, dims=1)
    maxs = maximum(reducedData, dims=1)
    return mean(maxs - mins)
end

# Z-score normalization function
function zscoreNormalize(dataMatrix)
    logdataMatrix = log.(dataMatrix)
    means = mean(logdataMatrix, dims=1)
    stds = std(logdataMatrix, dims=1)
    return (logdataMatrix .- means) ./ stds
end

# Min-Max scaling function
function minMaxScale(dataMatrix)
    mins = minimum(dataMatrix, dims=1)
    maxs = maximum(dataMatrix, dims=1)
    return (dataMatrix .- mins) ./ (maxs - mins)
end




# Normalize the data
normalizedData = zscoreNormalize(testpop)
normalizedData = minMaxScale(testpop)

# Reduce dimensions using PCA
reducedData = applyPCA(normalizedData)

# Calculate spread in reduced dimensions
spreadMetric = calculateSpread(reducedData)


popdf = make_pop_dataframe(test_population, allconstraints)

ranges = []
for constraint in allconstraints
    conrange = constraint.max - constraint.min
    push!(ranges, conrange)
end
    
normdf = popdf./ranges'

# using Plots
using StatsPlots
@df normdf boxplot(cols(), xticks = (1:activelength(allconstraints), propertynames(allconstraints)), outliers=false, legend = false, title = "Normalized Population Distributions", ylabel = "Normalized Value", xlabel = "Parameter", size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

#for each population size [1000, 5000, 10000, 15000, 20000, 25000, 30000], make the boxplot and add as subplot
poplots = []
for popsize in [1000, 5000, 10000, 20000, 50000, 100000, 200000]
    test_population = generate_population(allconstraints, popsize)
    popdf = make_pop_dataframe(test_population, allconstraints)
    normdf = popdf./ranges'
    pl = @df normdf boxplot(cols(), title = "$popsize Individuals",xticks = (1:activelength(allconstraints), propertynames(allconstraints)), outliers=false, legend = false, ylabel = "Normalized Value")
    push!(poplots, pl)
end
plot(poplots..., layout = (4, 2), size = (1200, 1200), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

using RDatasets
school = RDatasets.dataset("mlmRev","Hsb82")
# @df school density(:MAch, group = :Sx)
@df school density(:MAch, group = (:Sx, :Sector), legend = :topleft)

@df ga_df density(:per, :amp, legend = :topleft)



per_speed = []
for per in ga_df.per
    if per <= mean(ga_df.per)
        push!(per_speed, "Slow")
    else
        push!(per_speed, "Fast")
    end
end
ga_df.per_speed = per_speed

pl = @df ga_df density(:ka2, group = :per_speed, legend = false)
pl = density(ga_df.ka1, group = ga_df[:,:per_speed], legend = false)


plotvec = []
for parameter_name in propertynames(ga_df[:,Between(:ka1, :kcat7)])
    pl = density(ga_df[:,parameter_name], group = ga_df[:,:per_speed], legend = false, title=parameter_name, trim=true)
    push!(plotvec, pl)
end
plot(plotvec..., layout = (4, 3), size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)


@df ga_df andrewsplot(:per_speed, cols(7:15))

@df ga_df corrplot(cols(6:10))



BLAS.set_num_threads(18)
FFTW.set_num_threads(18)

@benchmark test4fixedGA(5000)
"""BenchmarkTools.Trial: 1 sample with 1 evaluation.
Single result which took 26.816 s (32.50% GC) to evaluate,
with a memory estimate of 92.33 GiB, over 547954601 allocations."""


#*FFTW testing 

ogprobjac = make_ODE_problem()

sol = solve_odeprob(ogprobjac, [6, 9, 10, 11, 12, 15, 16])



Amem_sol = map(sum, sol.u)

@benchmark getFrequencies(Amem_sol)
fft_array = getFrequencies(Amem_sol)
fft_array_downsampled = getFrequencies(Amem_sol[1:2:end])
idx_max, vals_max, idx_min, vals_min = OscTools.findextrema(fft_array; height = 1e-2, distance = 5)
idx_max_down, vals_max_down, idx_min_down, vals_min_down = OscTools.findextrema(fft_array_downsampled; height = 1e-2, distance = 5)

setdiff(vals_max, vals_max_down)

fft_array = @view Amem_sol[1:cld(length(Amem_sol),2)]

@benchmark OscTools.getFrequencies!($fft_array, $Amem_sol)

@benchmark FitnessFunction($Amem_sol, $sol.t)

FitnessFunction(Amem_sol, sol.t)


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