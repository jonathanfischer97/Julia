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
    gr()
    using Base.Threads
end

const PARAM_NAMES = ["ka1", "kb1", "kcat1", "ka2", "kb2", "ka3", "kb3", "ka4", "kb4", "ka7", "kb7", "kcat7", "DF"]

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 DF
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase, P = Phosphatase, A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*DF,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*DF,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

    #previously excluded reactions, all possible combinations possible in vitro
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*DF,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*DF,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*DF,kb3), LpA + LK <--> LpAKL
    (ka4*DF,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end  

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
    #? Parameter list
    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #? Initial condition list
    usym = [:L => 0.0, :Lp => 3.0, :K => 0.5, :P => 0.3, :A => 2.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    prob = ODEProblem(fullrn, u0, tspan, p)
    sol = solve(prob, saveat=0.1, save_idxs=1)

    #? Plot the results
    # plot(sol)
end




#! Helper functions for cost function ## 
begin
    """Get summed difference of peaks in the frequency domain"""
    function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) 
        idxarrLen = length(indexes)
        sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(idxarrLen-1))
        sum_diff += arrayData[indexes[end]]
        return sum_diff
    end

    """Get summed average standard deviation of peaks in the frequency domain"""
    function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64})#, window_ratio::Float64) #get average standard deviation of fft peak indexes
        arrLen = length(arrayData)
        window = 1 #max(1, round(Int, window_ratio * arrLen))
        sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in peakindxs)
        return sum_std / length(peakindxs)
    end

    """Return normalized FFT of solution vector"""
    function getFrequencies(y::Vector{Float64})
        res = abs.(rfft(y))
        return res ./ cld(length(y), 2) #normalize amplitudes
    end

    """Calculates the period and amplitude of each individual in the population"""
    function getPerAmp(sol::ODESolution)
        # Find peaks and calculate amplitudes and periods
        indx_max, vals_max = findmaxima(sol.u, 1)
        indx_min, vals_min = findminima(sol.u, 1)

        if length(indx_max) < 2 || length(indx_min) < 2
            return 0., 0.
        else
            # Calculate amplitudes and periods
            @inbounds pers = [sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1)]
            @inbounds amps = [(vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min))]

            # Calculate means of amplitudes and periods
            per = mean(pers)
            amp = mean(amps)

            return per, amp
        end
    end

    """Cost function to be plugged into eval_fitness wrapper"""
    function CostFunction(Y::ODESolution)
        #get the fft of the solution
        fftData = getFrequencies(Y.u)
        fftindexes = findmaxima(fftData,1)[1] #get the indexes of the peaks in the fft
        timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the fft
        if length(fftindexes) < 2 || length(timeindexes) < 2 #if there are no peaks, return 0
            return 0.0, 0.0, 0.0
        end
        std = getSTD(fftindexes, fftData) #get the standard deviation of the peaks
        diff = getDif(fftindexes, fftData) #get the difference between the peaks

        # Compute the period and amplitude
        period, amplitude = getPerAmp(Y)

        # Return cost, period, and amplitude as a tuple
        return -std + diff, period, amplitude
    end
end


"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(p,  prob::ODEProblem)
    Y = nothing
    try 
        Y = solve(remake(prob, p=p), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
        if Y.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(Y) for x in array) || any(x==1 for array in isless.(Y, 0.0) for x in array)
            return 1.0
        end
    catch e 
        if e isa DomainError #catch domain errors
            return 1.0
        else
            rethrow(e) #rethrow other errors
        end
    end
    fitness, period, amplitude = CostFunction(Y)
    return (fit = -fitness, per = period, amp = amplitude)
end


@bprofile eval_fitness_catcherrors(p, prob)





#! Helper functions for evaluating the contour of 2D parameter space ##
"""Find the indices of the parameters in the PARAM_NAMES array"""
function find_param_indices(param_combination)
    p1idx = findfirst(isequal(param_combination[1]), PARAM_NAMES)
    p2idx = findfirst(isequal(param_combination[2]), PARAM_NAMES)
    return p1idx, p2idx
end

"""Update the parameters inplace in the newparams array"""
function update_params!(newvals, newparams, p1idx, p2idx)
    newparams[p1idx] = newvals[1]
    newparams[p2idx] = newvals[2]
end



"""ParamRange struct for defining parameter ranges"""
struct ParamRange
    min::Float64
    max::Float64
end

"""Parameter constraint function, returns a dictionary of parameter names and their ranges"""
function define_parameter_constraints()
    # Define parameter constraint ranges
    ka_min, ka_max = -3.0, 1.0  # uM^-1s^-1, log scale
    kb_min, kb_max = -3.0, 3.0  # s^-1, log scale
    kcat_min, kcat_max = -3.0, 3.0 # s^-1, log scale
    df_min, df_max = 1.0, 5.0 # for DF, log scale

    param_values = OrderedDict(
        "ka1" => ParamRange(ka_min, ka_max),
        "kb1" => ParamRange(kb_min, kb_max),
        "kcat1" => ParamRange(kcat_min, kcat_max),
        "ka2" => ParamRange(ka_min, ka_max),
        "kb2" => ParamRange(kb_min, kb_max),
        "ka3" => ParamRange(ka_min, ka_max),
        "kb3" => ParamRange(kb_min, kb_max),
        "ka4" => ParamRange(ka_min, ka_max),
        "kb4" => ParamRange(kb_min, kb_max),
        "ka7" => ParamRange(ka_min, ka_max),
        "kb7" => ParamRange(kb_min, kb_max),
        "kcat7" => ParamRange(kcat_min, kcat_max),
        "DF" => ParamRange(df_min, df_max)
    )
    return param_values
end


#* Define the function to evaluate the 2D parameter space for oscillations
function evaluate_2D_parameter_space(paramrange_pair, prob::ODEProblem, paramrange_dict::OrderedDict=define_parameter_constraints(); steps=300)

    # Get indices of parameters
    p1idx, p2idx = find_param_indices(paramrange_pair)
    
    # Get evaluation function
    evalfunc(newparams) = eval_fitness_catcherrors(newparams, prob)
    
    # Get ranges of parameters
    p1_range = exp10.(range(paramrange_dict[paramrange_pair[1]].min, paramrange_dict[paramrange_pair[1]].max, length=steps))
    p2_range = exp10.(range(paramrange_dict[paramrange_pair[2]].min, paramrange_dict[paramrange_pair[2]].max, length=steps))

    # Get dimensions of parameter space
    p1_length = length(p1_range)
    p2_length = length(p2_range)

    # Create arrays to store results, oscillatory and non-oscillatory points
    oscillation_scores = Array{Float64}(undef, p1_length, p2_length)
    periods = Array{Float64}(undef, p1_length, p2_length)

    # Create progress bar
    innerprogress = Progress(p1_length*p2_length, dt = 0.1, desc="Evaluating $paramrange_pair parameter space... ", color=:red)

    # Evaluate parameter space
    @threads for i in 1:p1_length
        newparams = copy(prob.p)
        for j in 1:p2_length
            update_params!((p1_range[i], p2_range[j]), newparams, p1idx, p2idx)
            results = evalfunc(newparams)
            oscillation_scores[j, i] = results.fit
            periods[j, i] = results.per #< -0.5 ? results.per : 0.0
            next!(innerprogress)
        end
    end

    return (fit = -oscillation_scores, per = periods, ranges = (p1_range, p2_range), names = paramrange_pair)
end

function unit_labeler(parameter_name)
    if parameter_name in ["ka1", "ka2", "ka3", "ka4", "ka7"]
        return "$parameter_name (uM⁻¹ s⁻¹)"
    elseif parameter_name in ["kb1", "kb2", "kb3", "kb4", "kb7", "kcat1", "kcat7"]
        return "$parameter_name (s⁻¹)"
    elseif parameter_name in ["DF"]
        return "V/A (nm)"
    else
        return "unitless"
    end
end


function plot_oscillation_contour(result)
    p1_range = result.ranges[1]
    p2_range = result.ranges[2]

    oscplot = contour(p1_range, p2_range, result.fit, xscale=:log10, yscale=:log10, 
                        title="Oscillation Scores", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Oscillation Index", color=:vik, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)

    frequencies = ifelse.(result.per .== 0, 0, 1 ./ result.per)
    freqplot = contour(p1_range, p2_range, frequencies, xscale=:log10, yscale=:log10, 
                        title="Frequencies", xlabel=unit_labeler(result.names[1]), ylabel=unit_labeler(result.names[2]), colorbar_title="Frequency (Hz)", color=:plasma, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    plot(oscplot, freqplot, layout=(2,1), size=(1000, 800))
end


testresult = evaluate_2D_parameter_space(("ka1", "DF"), prob; steps =300)
testpers = testresult.per

#* Plot oscillatory regions of 2D parameter space with contour or heatplot
testplot = plot_oscillation_contour(testresult)
savefig(testplot, joinpath(@__DIR__, "testplot_ka2_DF.png"))




function evaluate_and_plot_all_2D_combinations(prob::ODEProblem, paramrange_dict::OrderedDict=define_parameter_constraints(); steps=300, savepath="./OscillatorPaper/FigureGenerationScripts/Figures/ContourPlots")
    # Get all combinations of parameters
    param_keys = collect(keys(paramrange_dict))
    param_combinations = combinations(param_keys, 2)
    
    # Create a dictionary to store all results
    result_dict = Dict()

    # Loop through each combination
    for combination in param_combinations
        # Evaluate the parameter space
        result = evaluate_2D_parameter_space(combination, prob, paramrange_dict; steps)

        # Store the result in the dictionary
        result_dict[combination[1]*"_"*combination[2]] = result

        # Plot the result
        comboplot = plot_oscillation_contour(result)

        #Save to a file
        savefig(comboplot, joinpath(savepath, "$(combination[1])_$(combination[2]).png"))
    end
    # Return the dictionary with all results
    return result_dict
end

result_dict = evaluate_and_plot_all_2D_combinations(prob; steps=300)


