begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/ODEProbMaker.jl")

    include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], (0.,1000.0), [])
new_u0 = ogprob.u0 .* 10
ogprob = remake(ogprob, u0 = new_u0)
# @benchmark solve($ogprob, saveat = 0.1, save_idxs = 1)

@benchmark solve($ogprob, saveat=0.1, save_idxs = 1)
ogsol = solve(ogprob, saveat=0.1, save_idxs = 1)
testfunc(ogprob) = solve(ogprob, saveat=0.1, save_idxs = 1)

using Cthulhu
using ProfileView
descend_code_warntype(testfunc, (ODEProblem,))
@code_warntype solve(ogprob, saveat=0.1, save_idxs = 1)
# plot(ogsol)
@code_warntype solve_for_fitness_peramp(ogprob)

@benchmark solve_for_fitness_peramp($ogprob)

@code_warntype CostFunction(ogsol)



#* Optimization of parameters to produce data for CSV
param_constraints = define_parameter_constraints(karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e3, 1e5))


# param_gaproblem = GAProblem(param_constraints,ogprob)

# param_record = run_GA(param_gaproblem)
# CSV.write("OscillatorPaper/FigureGenerationScripts/initialconditions.csv", param_record)




# """Compare the allowed solution space for when each triplet of parameters is fixed"""
# function reachability_analysis(constraints::ConstraintType, prob::ODEProblem) 
#     #* Get all combinations of fixed pairs
#     fixed_triple_combos = combinations(constraints.ranges, 3)

#     combos_length = length(collect(fixed_triple_combos))
#     combos_names = [string(triplet[1].name, "_", triplet[2].name, "_", triplet[3].name) for triplet in fixed_triple_combos]

#     loopprogress = Progress(length(combos_length), desc ="Looping thru fixed pairs: " , color=:blue)

#     #* Make a results DataFrame where fixedpair => (volume, avg_period, avg_amplitude, num_oscillatory_points)
#     results_df = DataFrame(FixedTriplet = combos_names, volume = Vector{Float64}(undef, combos_length), periodstats = Vector{Vector{Float64}}(undef, combos_length), amplitudestats = Vector{Vector{Float64}}(undef, combos_length), num_oscillatory_points = Vector{Int}(undef, combos_length))

#     for (i,fixedpair) in enumerate(fixed_triple_combos)
#         @info "Fixed input pair: $(fixedpair[1].name), $(fixedpair[2].name)"

#         #* Make a copy of the constraints to modify
#         variable_constraints = deepcopy(constraints)

#         #* Remove the fixed parameters from the variable constraints
#         filter!(x -> x.name != fixedpair[1].name && x.name != fixedpair[2].name, variable_constraints.ranges)

#         fixed_ga_problem = GAProblem(variable_constraints, prob)

#         fixedpair_idxs = find_indices(fixedpair, constraints.ranges) # Get the indices of the fixed inputs.

#         #* Create a closure for the fitness function that includes the fixed inputs
#         make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixedpair, fixedpair_idxs)

#         #* Run the optimization function to get the evaluated points
#         oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 8) #TODO: outputting the same number of points for multiple pairs

 
#         #* Calculate the number of oscillatory points
#         num_points = length(oscillatory_points_df.ind)
       

#         # #* Calculate the average period and amplitude
#         # periodstats = quantile(oscillatory_points_df.per, [0.0, 0.25, 0.5, 0.75, 1.0])
        

#         # amplitudestats = quantile(oscillatory_points_df.amp, [0.0, 0.25, 0.5, 0.75, 1.0])
        

#         #* Store the results for the given fixed parameter combination
#         results_df[i, 2:end] .= (volume, oscillatory_points_df.per, oscillatory_points_df.amp, num_points)

#         next!(loopprogress)
#     end
#     return results_df
# end



# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int})
    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if i in triplet_idxs
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_triplet[fixed_inputs_inserted + 1]
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[i - fixed_inputs_inserted]
            end
        end

        return evalfunc(new_input, prob)
    end
    return fitness_function
end



"""Defines logspace function for sampling parameters"""
logrange(start, stop, length) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))



#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; rangelength = 4) #TODO add progress bar
    variable_constraints = deepcopy(constraints)
    fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    # fixed_values1 = range(fixedtrip[1].min, stop = fixedtrip[1].max, length = 5)
    fixed_values1 = logrange(fixedtrip[1].min, fixedtrip[1].max, rangelength)
    # fixed_values2 = range(fixedtrip[2].min, stop = fixedtrip[2].max, length = 5)
    fixed_values2 = logrange(fixedtrip[2].min, fixedtrip[2].max, rangelength)
    # fixed_values3 = range(fixedtrip[3].min, stop = fixedtrip[3].max, length = 5)
    fixed_values3 = logrange(fixedtrip[3].min, fixedtrip[3].max, rangelength)

    results_df = DataFrame(param1 => Vector{Float64}(undef, rangelength^3), param2 => Vector{Float64}(undef, rangelength^3), param3 => Vector{Float64}(undef, rangelength^3), 
                        "average_period" => Vector{Float64}(undef, rangelength^3), "maximum_period"=>Vector{Float64}(undef, rangelength^3), "minimum_period"=>Vector{Float64}(undef, rangelength^3),
                        "average_amplitude" => Vector{Float64}(undef, rangelength^3), "maximum_amplitude"=>Vector{Float64}(undef, rangelength^3), "minimum_amplitude"=>Vector{Float64}(undef, rangelength^3))

    
    #* make folder to hold all the csv files 
    path = mkpath("OscillatorPaper/FigureGenerationScripts/3FixedCSVs/$(param1)_$(param2)_$(param3)")
    i = 1
    for val1 in fixed_values1
        for val2 in fixed_values2
            for val3 in fixed_values3
                fixed_values = [val1, val2, val3]
                @info fixed_values
                make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs)
                oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 5) 

                if isempty(oscillatory_points_df)
                    results_df[i, :] = (val1, val2, val3, NaN, NaN, NaN, NaN, NaN, NaN)
                    continue
                else
                    average_period::Float64 = mean(oscillatory_points_df.per)
                    maximum_period::Float64 = maximum(oscillatory_points_df.per)
                    minimum_period::Float64 = minimum(oscillatory_points_df.per)
    
                    average_amplitude::Float64 = mean(oscillatory_points_df.amp)
                    maximum_amplitude::Float64 = maximum(oscillatory_points_df.amp)
                    minimum_amplitude::Float64 = minimum(oscillatory_points_df.amp)
    
                    # push!(results_df, (val1, val2, val3, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude))
                    results_df[i, :] = (val1, val2, val3, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude)
                    i += 1
                
                    #* insert the fixed params into each ind of oscillatory_points_df
                    for (rowid,row) in enumerate(eachrow(oscillatory_points_df))
                        for (j,fixedidx) in enumerate(fixedtrip_idxs)
                            if fixedidx <= length(row.ind)
                                insert!(oscillatory_points_df[rowid,:ind], fixedtrip_idxs[j], fixed_values[j])
                            else
                                push!(oscillatory_points_df[rowid,:ind], fixed_values[j])
                            end
                        end
                    end
                    @assert length(oscillatory_points_df.ind[1]) == 13
                    return oscillatory_points_df
                    #* split parameter values into separate columns and add initial conditions
                    split_dataframe!(oscillatory_points_df, prob)
                    CSV.write(path*"/$(round(val1; digits = 2))_$(round(val2;digits = 2))_$(round(val3; digits=2)).csv", oscillatory_points_df)
                end
            end
        end
    end
    return results_df
end

param_triplet = ["kcat1", "kcat7", "DF"]

results_df = fixed_triplet_csv_maker(param_triplet..., param_constraints, ogprob)
testdf = deepcopy(results_df)

split_dataframe!(testdf, ogprob)

CSV.write("OscillatorPaper/FigureGenerationScripts/fixed_triplet_results-$(param_triplet[1]*param_triplet[2]*param_triplet[3]).csv", results_df)

for rowid in eachrow(testdf)
    println(rowid)
end

@code_warntype fixed_triplet_csv_maker("ka1", "ka2", "ka3", param_constraints, ogprob)




#! TESTING GA FUNCTIONALITY
test_gaproblem = GAProblem(param_constraints, ogprob)

ga_results = run_GA(test_gaproblem; population_size = 1000, iterations = 5)

split_dataframe!(ga_results, ogprob)

test_results

@time run_GA(test_gaproblem; population_size = 1000, iterations = 5)
sort!(test_results, :fit, rev=false)


plotsol(row) = plotsol(row, test_results, ogprob)

plotsol(1)


test_fitness(row) = eval_param_fitness(test_results.ind[row], ogprob)
test_fitness(1)

reogprob = remake(ogprob, p=test_results.ind[1])
testsol = solve(reogprob, saveat = 0.01, save_idxs = 1)
CostFunctionTest(testsol)

plot(testsol)
CostFunction(testsol)
@benchmark CostFunction($testsol)


 
#Test the fitness function with scaled amplitudes
testsolx10 = testsol.u .* 10

CostFunction(testsolx10)

"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunctionTest(sol::ODESolution)::Vector{Float64}
    # tstart = 50 #iterations, = 5 seconds #TODO look if last half of sol is constant, if so, cut it off
    # trimsol = sol[tstart:end] #* get the solution from the clean start time to the end
    normsol = normalize_time_series(sol[1,:]) #* normalize the solution
    #*get the fft of the solution
    fftData = getFrequencies(normsol)
    fft_peakindexes, fft_peakvals = findmaxima(fftData,5) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(sol[1,:],5) #* get the times of the peaks in the fft
    @debug "FFT peak values: $fft_peakvals"
    @debug "Time Peak indexes: $time_peakindexes"
    @debug "Time Peak values: $time_peakvals"
    if length(fft_peakindexes) < 1 #* if there are no peaks in either domain, return 0
        @debug "No peaks found in FFT"
        return [0.0, 0.0, 0.0]
    elseif length(time_peakindexes) < 3
        @debug "Not enough peaks found in time domain"
        return [0.0, 0.0, 0.0]
    elseif maximum(peakproms(time_peakindexes,sol[1,:])[2]) < 0.1
        @debug "No peaks found in time domain with prominence > 0.1"
        return [0.0, 0.0, 0.0]
    else
        std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
        @debug "std = $std"
        diff = getDif(fft_peakvals) #* get the summed difference between the peaks in frequency domain
        @debug "diff = $diff"

        #* Compute the period and amplitude
        period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)

        return [-std - diff, period, amplitude]
    end

    #* if cost is too high, return 0, no oscillations
    # if -std-diff > -0.1
    #     return [0.0, 0.0, 0.0]
    # else 
    #     period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)
    #     return [-std - diff, period, amplitude]
    # end
end

function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64})
    # Find peaks of the minima too
    indx_min, vals_min = findminima(sol[1,:], 1)

    if length(indx_max) < 2 || length(indx_min) < 2
        return 0.0, 0.0
    end

    # Calculate amplitudes
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))

    # Calculate periods based on the maxima
    @inbounds pers_max = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    
    # Normalize periods by total time span
    total_time_span = sol.t[end] - sol.t[1]
    normalized_pers_max = [p / total_time_span for p in pers_max]

    return mean(normalized_pers_max), mean(amps)
end


getFrequencies(testsol.u)
getFrequencies(testsolx10)


normsol = normalize_time_series(testsol[1,:]) #* normalize the solution
@benchmark normalize_time_series!($testsol[1,:])
normsol = normalize(testsol[1,:], 1) #* normalize the solution
#*get the fft of the solution
fftData = getFrequencies(normsol)
fft_peakindexes, fft_peakvals = findmaxima(fftData,5) #* get the indexes of the peaks in the fft
time_peakindexes, time_peakvals = findmaxima(testsol[1,:],5) #* get the times of the peaks in the fft
# timeproms = peakproms(time_peakindexes, normsol)[2]
# std(timeproms)
if length(fft_peakindexes) < 3 || length(time_peakindexes) < 3 #* if there are no peaks in either domain, return 0
    return [0.0, 0.0, 0.0]
end
standev = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
peakdiff = getDif(fft_peakvals) #* get the summed difference between the peaks in frequency domain

#* Compute the period and amplitude
period, amplitude = getPerAmp(testsol, time_peakindexes, time_peakvals)


@benchmark getFrequencies(ogsol[1,:])

plot(fftData)
plot(testsol)

#< Bugs to fix
#* 1. The fitness function is not working properly. Optimize, reduce fine tuning, and test it for edge cases against expected results
#! * 2. Fix the "FFTW can't make plan" error 
#! * 3. Logscale projection isn't working, fix it. 
#! 4. Save the optimized parameters, not just the evaluation values 
#* 5. Run GA through debugger to see the sequence of selection, recombination
#! 6. Print out the initial conditions in the CSV 
#* 7. Make test suite for the fitness function. Orthogonal tests will run through all optimized solutions and classify them as correct, false negative, false positive.
#* 8. Save all variables for selected solutions. IDK Maggie mentioned it 



newp = [0.001741312,	0.002189009,	218.8342196,	4.333715377,	383.0505446,	
0.858290292,	19.19847775,	1.249452352,	0.045457933,	0.330913943,	0.007674134,	0.006808512,	75228.30051]

newprob = remake(ogprob, p = newp)

newsol = solve(newprob, saveat = 0.01, save_idxs = 1)
plot(newsol)

CostFunction(newsol)


testarray = rand(10000)
testfft = getFrequencies(testarray)
@benchmark getFrequencies($testarray)


# load CSV into DataFrame 
testdf = CSV.read("/Users/jonathanfischer/Desktop/PhD_ThesisWork/Julia/OscillatorPaper/FigureGenerationScripts/fixed_triplet_results-kcat1kcat7DF.csv", DataFrame)

# exclude rows with extreme values 
testdf = testdf[.!isnan.(testdf.average_period), :]

plot(testdf.kcat1, testdf.kcat7, testdf.DF, st = :scatter, xlabel = "kcat1", ylabel = "kcat7", zlabel="DF", title = "Average Period vs. Average Amplitude", legend = false)





# Clean the data by removing extreme values
df = filter(row -> !any(isinf.(row)) && !any(isnan.(row)), testdf)

df = testdf
# Convert to log scale
log_kcat1 = log10.(df.kcat1)
log_kcat7 = log10.(df.kcat7)
log_DF = log10.(df.DF)
sizes = df.average_amplitude .+ 1 # Add 1 to avoid zero size
colors = df.average_period

# Plotting function
function create_3d_scatter_with_shadows(angle, x, y, z, sizes, colors)
    p = scatter3d(x, y, z, markersize=sizes, color=cgrad(:viridis, colors, rev=true), legend=false, alpha=0.4, markerstrokewidth=0)
    display(p)
    # Project the points onto each plane and connect them with lines
    for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, colors)
        scatter3d!([xi], [yi], [z[begin]], markersize=[si], color=cgrad(:viridis, [ci], rev=true), alpha=0.4)
        scatter3d!([xi], [y[begin]], [zi], markersize=[si], color=cgrad(:viridis, [ci], rev=true), alpha=0.4)
        scatter3d!([x[begin]], [yi], [zi], markersize=[si], color=cgrad(:viridis, [ci], rev=true), alpha=0.4)
    end

    # Set axis labels
    xlabel!("log10(kcat1)")
    ylabel!("log10(kcat7)")
    zlabel!("log10(DF)")

    # Set view angle
    plot!(camera=(10, angle))

    # Set title
    title!("Log-Scaled Parameters with Projections & Shadows (Angle: ) - Size & Color: Average Amplitude & Period")

    # Display the plot
    # display(current())
    
end

# Plot from different angles with log scaling and projections
angles = [0, 30, 60, 90, 120, 150]
for ang in angles
    create_3d_scatter_with_shadows(ang, log_kcat1, log_kcat7, log_DF, sizes, colors)
end
