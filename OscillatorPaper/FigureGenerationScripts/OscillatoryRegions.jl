begin 
    using Plots 
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    # using Unitful
    # using Unitful: µM, nm, s
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
    # plotlyjs()
    # gr()
    using Base.Threads
end

# import the Catalyst model "fullrn"
include("../../UTILITIES/ReactionNetwork.jl")

# import the cost function and other evaluation functions
include("../../UTILITIES/EvaluationFunctions.jl")

# import the genetic algorithm and associated functions
include("../../UTILITIES/GA_functions.jl")

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
    #? Parameter list
    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #? Initial condition list
    usym = [:L => 3.0, :Lp => 0.0, :K => 0.5, :P => 0.3, :A => 2.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    fullprob = ODEProblem(fullrn, u0, tspan, p)
    sol = solve(fullprob, saveat=0.1, save_idxs=1)

    #? Plot the results
    plot(sol)
end





#! Helper functions for evaluating the contour of 2D parameter space ##
"""Update the parameters inplace in the newparams array"""
function update_params!(newvals, newparams::Vector{Float64}, p1idx, p2idx)
    newparams[p1idx] = newvals[1]
    newparams[p2idx] = newvals[2]
end




function get_range(paramrange::Symbol, constraints::ParameterConstraints, length::Int)
    return exp10.(range(getproperty(constraints,paramrange).min, getproperty(constraints,paramrange).max, length=length))
end

function get_range(ICrange::Symbol, constraints::InitialConditionConstraints, length::Int)
    return range(getproperty(constraints,ICrange).min, getproperty(constraints,ICrange).max, length=length)
end



#! Define the function to evaluate the 2D solution space for oscillations
function evaluate_2D_solution_space(paramrange_pair::Vector{Symbol}, prob::ODEProblem, constraints::ConstraintType; steps=300)
    # Get string names from pair
    p1name = getproperty(constraints, paramrange_pair[1]).name
    p2name = getproperty(constraints, paramrange_pair[2]).name

    constraints isa ParameterConstraints ? NAMES = PARAM_NAMES : NAMES = VAR_NAMES #TODO Make this generic

    # Get indices of parameters
    p1idx, p2idx = find_indices(paramrange_pair, NAMES)
    
    # Get evaluation function
    #TODO Make this generic, too much type checking
    evalfunc = constraints isa ParameterConstraints ? eval_param_fitness : eval_ic_fitness

    
    # Get ranges of parameters
    p1_range = get_range(paramrange_pair[1], constraints, steps)
    p2_range = get_range(paramrange_pair[2], constraints, steps)

    # Get dimensions of parameter space
    p1_length = length(p1_range)
    p2_length = length(p2_range)

    # Create arrays to store results, oscillatory and non-oscillatory points
    oscillation_scores = Array{Float64}(undef, p1_length, p2_length)
    periods = Array{Float64}(undef, p1_length, p2_length)

    # Create progress bar
    innerprogress = Progress(p1_length*p2_length, dt = 0.1, desc="Evaluating $p1name : $p2name solution space... ", color=:red)


    @threads for i in 1:p1_length
        newparams = copy(prob.u0) #TODO Make this generic
        for j in 1:p2_length
            update_params!((p1_range[i], p2_range[j]), newparams, p1idx, p2idx)
            results = evalfunc(newparams,prob)
            oscillation_scores[j, i] = results[1]
            periods[j, i] = results[1] < -0.2 ? results[2] : 0.0
            next!(innerprogress)
        end
    end
    return (fit = -oscillation_scores, per = periods, ranges = (p1_range, p2_range), names = (p1name, p2name))
end




function unit_labeler(name)
    if name in ["ka1", "ka2", "ka3", "ka4", "ka7"]
        return "$name (μM⁻¹ s⁻¹)"
    elseif name in ["kb1", "kb2", "kb3", "kb4", "kb7", "kcat1", "kcat7"]
        return "$name (s⁻¹)"
    elseif name in ["DF"]
        return "V/A (nm)"
    else #is IC concentration
        return "$name (μM)"
    end
end


function plot_oscillation_contour(result; xscaling=:log10, yscaling=:log10)
    var1_range = result.ranges[1]
    var2_range = result.ranges[2]

    oscplot = contour(var1_range, var2_range, result.fit, xscale=xscaling, yscale=yscaling, 
                        title="Oscillation Scores", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Oscillation Index", color=:vik, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)

    frequencies = ifelse.(result.per .== 0, 0, 1 ./ result.per)
    freqplot = contour(var1_range, var2_range, frequencies, xscale=xscaling, yscale=yscaling, 
                        title="Frequencies", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Frequency (Hz)", color=:plasma, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    plot(oscplot, freqplot, layout=(2,1), size=(1000, 1000))
end

function plot_oscillation_heatmap(result; xscaling=:log10, yscaling=:log10)
    var1_range = collect(result.ranges[1])
    var2_range = collect(result.ranges[2])

    oscplot = heatmap(var1_range, var2_range, result.fit, xscale=xscaling, yscale=yscaling,
                        title="Oscillation Scores", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Oscillation Index", color=:vik, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)

    frequencies = ifelse.(result.per .== 0, 0, 1 ./ result.per)
    freqplot = heatmap(var1_range, var2_range, frequencies, xscale=xscaling, yscale=yscaling,
                        title="Frequencies", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Frequency (Hz)", color=:plasma, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    plot(oscplot, freqplot, layout=(2,1), size=(1000, 1000))
end

function plot_oscillation_surface(result; xscaling=:log10, yscaling=:log10)
    plotlyjs()
    var1_range = collect(result.ranges[1])
    var2_range = collect(result.ranges[2])

    oscplot = surface(var1_range, var2_range, result.fit, xscale=xscaling, yscale=yscaling,
                        title="Oscillation Scores", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Oscillation Index", color=:vik, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)

    frequencies = ifelse.(result.per .== 0, 0, 1 ./ result.per)
    freqplot = surface(var1_range, var2_range, frequencies, xscale=xscaling, yscale=yscaling,
                        title="Frequencies", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Frequency (Hz)", color=:plasma, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    plot(oscplot, freqplot, layout=(2,1), size=(1000, 1000))
    gr()
end


function plot_oscillation_scatter(result; xscaling=:log10, yscaling=:log10, step = 5)
    var1_range = result.ranges[1]
    sliced_var1_range = var1_range[1:step:end]

    var2_range = result.ranges[2]
    sliced_var2_range = var2_range[1:step:end]

    fitness_array = result.fit
    sliced_fitness_array = fitness_array[1:step:end,1:step:end]

    frequency_array = ifelse.(result.per .== 0, 0, 1 ./ result.per)
    sliced_frequency_array = frequency_array[1:step:end,1:step:end]


    # Initialize the plot with the first column or row
    oscplot = scatter(xlims=(0,var1_range[end]),ylims=(0,var2_range[end]), zlims = (minimum(sliced_fitness_array), maximum(sliced_fitness_array)), xscale=xscaling, yscale=yscaling,
                        title="Oscillation Scores", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Oscillation Index", label=false)

    # Do the same for the frequencies plot
    freqplot = scatter(xlims=(0,var1_range[end]),ylims=(0,var2_range[end]), zlims = (minimum(sliced_frequency_array), maximum(sliced_frequency_array)), xscale=xscaling, yscale=yscaling,
                        title="Frequencies", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Frequency (Hz)")


    # # Add the rest of the columns or rows iteratively
    # for j in 1+step:step:length(var2_range)
    #     scatter!(oscplot, sliced_var1_range, fill(var2_range[j],length(sliced_var1_range)), marker_z = fitness_array[j,1:step:end], color = :vik, label=false)
    # end

    for (i,var1) in enumerate(sliced_var1_range)
        for (j,var2) in enumerate(sliced_var2_range)
            if sliced_fitness_array[j,i] < 0.2
                scatter!(oscplot, [var1], [var2], color = :white, label=false, markeralpha = 0.0)
                scatter!(freqplot, [var1], [var2], color = :white, label=false, markeralpha = 0.0)
            else
                scatter!(oscplot, [var1], [var2], marker_z = [sliced_fitness_array[j,i]], color = :amp, label=false)
                scatter!(freqplot, [var1], [var2], marker_z = [sliced_frequency_array[j,i]], color = :lighttest, label=false)
            end
        end
    end

    # Combine both plots
    plot(oscplot, freqplot, layout=(2,1), size=(1000, 1000), bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
end



icranges = define_initialcondition_constraints()

testresult = evaluate_2D_solution_space([:L, :A], fullprob, icranges; steps =300)

#* Plot oscillatory regions of 2D parameter space with contour or heatplot
testcontour = plot_oscillation_contour(testresult; xscaling=:identity, yscaling=:identity)
testscatter = plot_oscillation_heatmap(testresult; xscaling=:identity, yscaling=:identity)
testsurface = plot_oscillation_surface(testresult; xscaling=:identity, yscaling=:identity)
testscatter = plot_oscillation_scatter(testresult; xscaling=:identity, yscaling=:identity)
savefig(testscatter, joinpath(@__DIR__, "testscatter_PIP_AP2.png"))




function evaluate_and_plot_all_2D_combinations(prob::ODEProblem, constraints::ConstraintType; steps=300, savepath="./OscillatorPaper/FigureGenerationScripts/Figures/ContourPlots")
    # Get all combinations of parameters
    param_keys = collect(keys(constraints.data))
    param_combinations = combinations(param_keys, 2)
    
    # Create a dictionary to store all results
    result_dict = Dict()

    # Loop through each combination
    for combination in param_combinations
        # Evaluate the parameter space
        result = evaluate_2D_solution_space(combination, prob, constraints; steps)

        # Get string names from pair
        name1, name2 = result.names

        # Store the result in the dictionary
        result_dict[name1*"_"*name2] = result

        # Plot the result depending on the type of constraints
        constraints isa ParameterConstraints ? comboplot = plot_oscillation_contour(result) : comboplot = plot_oscillation_contour(result; xscaling=:identity, yscaling=:identity)

        #Save to a file
        savefig(comboplot, joinpath(savepath, "$(name1)_$(name2).png"))
    end
    # Return the dictionary with all results
    return result_dict
end

icranges = define_initialcondition_constraints(;lipidrange=(0.1, 20.0), kinaserange=(0.1,5.0), phosphataserange=(0.1,5.0), ap2range=(0.1,20.0))
IC_result_dict = evaluate_and_plot_all_2D_combinations(fullprob, icranges; steps=1000)

