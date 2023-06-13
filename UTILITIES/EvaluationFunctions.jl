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
module EvaluationFunctions #< MODULE START

using DifferentialEquations: ODEProblem, ODESolution, solve #* for ODESolution type
using FFTW: rfft #* for FFT

#* Exported functions #####
export getPerAmp, CostFunction
#*#######



"""Get summed average difference of peaks in the frequency domain"""
function getDifAvg(peakvals::Vector{Float64}) #todo: compressed the range of values
    return (peakvals[1] - peakvals[end]) / length(peakvals)
end

#< START 
"""Get summed difference of peaks in the frequency domain"""
function getDif(peakvals::Vector{Float64})
    return peakvals[1] - peakvals[end]
end

"""Get summed average standard deviation of peaks in the frequency domain"""
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window::Int = 1)#, window_ratio::Float64) #get average standard deviation of fft peak indexes
    arrLen = length(fft_arrayData)
    sum_std = @inbounds sum(std(fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations
    return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end

"""Return normalized FFT of solution vector"""
function getFrequencies(timeseries::Vector{Float64}) #todo fix normalization or something 
    res = abs.(rfft(timeseries))
    return res ./ cld(length(timeseries), 2) #* normalize amplitudes
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64})
    #* Find peaks of the minima too 
    indx_min, vals_min = findminima(sol[1,:], 1)

    if length(indx_max) < 1 || length(indx_min) < 1 #todo need to fix this, either keep or move check into cost function
        return 0.0, 0.0
    end

    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))
    # @inbounds amps = 0.5 .* (vals_max .- vals_min)[1:min(length(indx_max), length(indx_min))]

    return mean(pers), mean(amps)
end

function getPerAmp(sol::ODESolution)
    #* Find peaks of the minima too 
    indx_max, vals_max = findmaxima(sol[1,:], 5)
    indx_min, vals_min = findminima(sol[1,:], 5)

    if length(indx_max) < 2 || length(indx_min) < 2
        return 0.0, 0.0
    end
    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))
    # @inbounds amps = 0.5 .* (vals_max .- vals_min)[1:min(length(indx_max), length(indx_min))]

    return mean(pers), mean(amps)
end

"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution)::Vector{Float64}
    tstart = 50 #iterations, = 5 seconds
    trimsol = sol[tstart:end] #* get the solution from the clean start time to the end
    # viewsol = sol.u #* get the view of the solution
    #*get the fft of the solution
    fftData = getFrequencies(trimsol[1,:])
    fft_peakindexes, fft_peakvals = findmaxima(fftData,1) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(trimsol[1,:],1) #* get the times of the peaks in the fft
    if length(fft_peakindexes) < 2 || length(time_peakindexes) < 2 #* if there are no peaks in either domain, return 0
        return [0.0, 0.0, 0.0]
    end
    std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
    diff = getDif(fft_peakvals) #* get the summed difference between the peaks in frequency domain

    #* Compute the period and amplitude
    period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)

    #* Return cost, period, and amplitude as a vector
    return [-std - diff, period, amplitude]
end


#! EVALUATION FUNCTIONS ## 
"""Evaluate the fitness of an individual with new parameters"""
function eval_param_fitness(params::Vector{Float64},  prob::ODEProblem)
    # remake with new parameters
    new_prob = remake(prob, p=params)
    return solve_for_fitness_peramp(new_prob)
end

"""Evaluate the fitness of an individual with new initial conditions"""
function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::ODEProblem)
    # remake with new initial conditions
    new_prob = remake(prob, u0=[initial_conditions; zeros(length(prob.u0)-length(initial_conditions))])
    return solve_for_fitness_peramp(new_prob)
end

"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function solve_for_fitness_peramp(prob)

    sol = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)

    return CostFunction(sol)
end


# """Utility function to solve the ODE and return the fitness"""
# function solve_for_fitness(prob::ODEProblem)

#     sol = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)

#     fitness = CostFunction(sol)[1]
#     return fitness
# end

end; #>MODULE END
