#! Helper functions for cost function ## 

# """Get summed difference of peaks in the frequency domain"""
# function getDif(peakvals::Vector{Float64}) #todo fix normalization
#     idxarrLen = length(peakvals)
#     sum_diff = @inbounds sum(peakvals[i] - peakvals[i+1] for i in 1:(idxarrLen-1))
#     # @info "Test sum diff of just first and last elements: $(peakvals[1] - peakvals[end])"
#     sum_diff += peakvals[end]
#     # @info "Old sum diff: $sum_diff"
#     return sum_diff
# end

# function getDif_bidirectional(peakvals::Vector{Float64})
#     idxarrLen = length(peakvals)

#     #* iterate from both ends simultaneously to deal with symmetry
#     sum_diff = @inbounds sum((peakvals[i] - peakvals[idxarrLen + 1 - i]) for i in 1:(idxarrLen ÷ 2)) 

#     return 2 * sum_diff #* multiply by 2 to account for the fact that we're only summing half of the differences
# end

"""Get summed average difference of peaks in the frequency domain"""
function getDifAvg(peakvals::Vector{Float64})
    return (peakvals[1] - peakvals[end]) / length(peakvals)
end

#< START 
"""Get summed difference of peaks in the frequency domain"""
function getDif(peakvals::Vector{Float64})
    return peakvals[1] - peakvals[end]
end

"""Get summed average standard deviation of peaks in the frequency domain"""
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window::Int = 2)#, window_ratio::Float64) #get average standard deviation of fft peak indexes
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
    indx_min, vals_min = findminima(sol.u, 5)

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

    if length(indx_max) < 3 || length(indx_min) < 3
        return 0.0, 0.0
    end
    #* Calculate amplitudes and periods
    @inbounds pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    @inbounds amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))
    # @inbounds amps = 0.5 .* (vals_max .- vals_min)[1:min(length(indx_max), length(indx_min))]

    return mean(pers), mean(amps)
end

"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution)::Tuple{Float64, Float64, Float64}
    #*get the fft of the solution
    fftData = getFrequencies(sol.u)
    fft_peakindexes, fft_peakvals = findmaxima(fftData,10) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(sol.u,5) #* get the times of the peaks in the fft
    if length(fft_peakindexes) < 2 || length(time_peakindexes) < 2 #* if there are no peaks in either domain, return 0
        return [0.0, 0.0, 0.0]
    end
    std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
    diff = getDif(fft_peakvals) #* get the summed difference between the peaks in frequency domain

    #* Compute the period and amplitude
    period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)

    #* Return cost, period, and amplitude as a vector
    return (-std - diff, period, amplitude)
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

# """Utility function to solve the ODE and return the fitness"""
# function solve_for_fitness(prob::ODEProblem)

#     sol = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)

#     fitness = CostFunction(sol)[1]
#     return fitness
# end

"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function solve_for_fitness_peramp(prob)

    sol = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)

    return CostFunction(sol)
end


# """Custom data structure to store the period and amplitude of each individual"""
# struct PeriodAmplitudes
#     peramps::Dict{Vector{Float64}, Tuple{Float64, Float64}}

#     PeriodAmplitudes() = new(Dict{Vector{Float64}, Tuple{Float64, Float64}}())
# end    

# """Helper function to update the period and amplitude of an individual, using in place in CostFunction"""
# function update_peramp!(tracker::PeriodAmplitudes, parameters::Vector{Float64}, values::Tuple{Float64, Float64})
#     tracker.peramps[parameters] = values
# end


# """Evaluate the fitness of an individual with new parameters and track periods and amplitudes"""
# function eval_param_fitness(params::Vector{Float64},  prob::ODEProblem, tracker::PeriodAmplitudes)
#     # remake with new parameters
#     new_prob = remake(prob, p=params)
#     return eval_fitness_catcherrors!(new_prob, tracker)
# end

# """Evaluate the fitness of an individual with new initial conditions and track periods and amplitudes"""
# function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::ODEProblem, tracker::PeriodAmplitudes)
#     # remake with new initial conditions
#     new_prob = remake(prob, u0=initial_conditions)
#     return eval_fitness_catcherrors!(new_prob, tracker)
# end

# """Cost function that also updates the period and amplitude in the tracker"""
# function eval_fitness_catcherrors!(prob::ODEProblem, peramp_tracker::PeriodAmplitudes)
#     Y = nothing
#     try 
#         Y = solve(prob, saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
#         if Y.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(Y) for x in array) || any(x==1 for array in isless.(Y, 0.0) for x in array)
#             return 1.0
#         end
#     catch e 
#         if e isa DomainError #catch domain errors
#             return 1.0
#         else
#             rethrow(e) #rethrow other errors
#         end
#     end
#     fitness, period, amplitude = CostFunction(Y)

#     # Update the additional values stored in the tracker
#     update_peramp!(peramp_tracker, p, (period, amplitude))

#     return -fitness
# end



