# Helper functions for cost function ## 

# """Get summed difference of peaks in the frequency domain"""
# function getDif(peakvals::Vector{Float64}) # fix normalization
#     idxarrLen = length(peakvals)
#     sum_diff = @inbounds sum(peakvals[i] - peakvals[i+1] for i in 1:(idxarrLen-1))
#     # @info "Test sum diff of just first and last elements: $(peakvals[1] - peakvals[end])"
#     sum_diff += peakvals[end]
#     # @info "Old sum diff: $sum_diff"
#     return sum_diff
# end

# function getDif_bidirectional(peakvals::Vector{Float64})
#     idxarrLen = length(peakvals)

#     # iterate from both ends simultaneously to deal with symmetry
#     sum_diff = @inbounds sum((peakvals[i] - peakvals[idxarrLen + 1 - i]) for i in 1:(idxarrLen ÷ 2)) 

#     return 2 * sum_diff # multiply by 2 to account for the fact that we're only summing half of the differences
# end


"""
# Module holding all evaluation functions for assesing oscillatory solutions
"""
# module EvaluationFunctions #< MODULE START

# using DifferentialEquations: ODEProblem, ODESolution, solve, remake #* for ODESolution type
# using FFTW: rfft #* for FFT

# #* Exported functions #####
# export getPerAmp, CostFunction, eval_ic_fitness, eval_param_fitness
# #*#######



# """Get summed average difference of peaks in the frequency domain"""
# function getDifAvg(peakvals::Vector{Float64}) #todo: compressed the range of values
#     return (peakvals[1] - peakvals[end]) / length(peakvals)
# end

#< START 
"""Get summed difference of peaks in the frequency domain"""
function getDif(peakvals::Vector{Float64})
    length(peakvals) == 1 ? peakvals[1] : peakvals[1] - peakvals[end]
end

"""Get summed average standard deviation of peaks in the frequency domain"""
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window_ratio::Int=100) #get average standard deviation of fft peak indexes
    arrLen = length(fft_arrayData)
    window = max(1,cld(arrLen,window_ratio)) #* window size is 1% of array length, or 1 if array length is less than 100
    sum_std = @inbounds sum(std(fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations
    return sum_std
    # return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end 

"""Return normalized FFT of solution vector. Modifies the solution vector in place"""
function getFrequencies(timeseries::Vector{Float64}) 
    res = abs.(rfft(timeseries))
    return res ./ cld(length(timeseries), 2) #* normalize by length of timeseries
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution)
    solu = sol[1,:] #* get the solution array out of ODESolution type

    indx_max, vals_max = findmaxima(solu, 5) #* find the peaks of the maxima with window size 5
    indx_min, vals_min = findminima(solu, 5)
    
    if length(indx_max) < 2 || length(indx_min) < 2
        return 0.0, 0.0
    end

    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    # @inbounds amps = ((solu[i] - solu[j])/2 for (i,j) in zip(indx_max,indx_min))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))

    return mean(pers), mean(amps)
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64})
    #* Find peaks of the minima too 
    indx_min, vals_min = findminima(sol[1,:], 5)

    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))

    return mean(pers), mean(amps)
end


"""Normalize a time series to have mean 0 and amplitude 1 before FFT"""
function normalize_time_series!(ts::Vector{Float64})::Vector{Float64}
    mu = mean(ts)
    amplitude = maximum(ts) - minimum(ts)
    ts .= (ts .- mu) ./ amplitude
end


"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution)::Vector{Float64}

    #* Check if last half of the solution array is steady state
    lasthalfsol = sol[cld(length(sol.t),2):end]
    if std(lasthalfsol[1,:]) < 0.01 
        return [0.0, 0.0, 0.0]
    end 

    #* Trim first 10% of the solution array to avoid initial spikes
    # tstart = cld(length(sol.t),10) 
    # trimsol = sol[tstart:end] 
    # trimsol = sol

    #* Get the solution array out of ODESolution type
    solarray = sol[1,:] 
    # time_peakindexes, time_peakvals = findmaxima(solarray,5) #* get the times of the peaks in the fft
    # time_peakproms = peakproms(time_peakindexes, solarray; minprom = 0.1)[1] #* get the peak prominences (amplitude of peak above surrounding valleys
    # if length(time_peakproms) < 2 #* if there are less than 2 prominent peaks in the time domain, return 0
    #     return [0.0, 0.0, 0.0]
    # end

    #* Normalize the solution array. WARNING: solarray is modified after this line
    normsol = normalize_time_series!(solarray)

    #* Get the rfft of the solution
    fftData = getFrequencies(normsol)
    # fft_peakindexes, fft_peakvals = findmaxima(fftData,5) #* get the indexes of the peaks in the fft
    fft_peakindexes, peakprops = findpeaks1d(fftData; height = 1e-2) #* get the indexes of the peaks in the fft
    fft_peakvals = peakprops["peak_heights"]
    if isempty(fft_peakindexes) #* if there is no signal in the frequency domain, return 0.0s
        return [0.0, 0.0, 0.0]
    else
        standard_deviation = getSTD(fft_peakindexes, fftData) #* get the summed standard deviation of the peaks in frequency domain
        sum_diff = getDif(fft_peakvals) #* get the summed difference between the first and last peaks in frequency domain
    
        #* Compute the period and amplitude
        period, amplitude = getPerAmp(sol)
    
        return [-standard_deviation - sum_diff, period, amplitude]
    end
end


#! EVALUATION FUNCTIONS ## 
"""Evaluate the fitness of an individual with new parameters"""
function eval_param_fitness(params::Vector{Float64},  prob::ODEProblem)
    #* remake with new parameters
    new_prob = remake(prob, p=params)
    return solve_for_fitness_peramp(new_prob)
end

"""Evaluate the fitness of an individual with new initial conditions"""
function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::ODEProblem)
    #* remake with new initial conditions
    new_prob = remake(prob, u0=[initial_conditions; zeros(length(prob.u0)-length(initial_conditions))])
    return solve_for_fitness_peramp(new_prob)
end

"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function solve_for_fitness_peramp(prob::ODEProblem)

    sol = solve(prob,saveat=0.1, save_idxs=1, verbose=false)
    # return CostFunction(sol)
    
    if sol.retcode == ReturnCode.Success
        return CostFunction(sol)
    else
        return [0.0, 0.0, 0.0]
    end
end

#>MODULE END
