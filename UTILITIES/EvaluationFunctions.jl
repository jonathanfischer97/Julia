#! Helper functions for cost function ## 

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

"""Cost function that also updates the period and amplitude in the tracker"""
function eval_fitness_catcherrors!(tracker::PeriodAmplitudes, p::Vector{Float64},  prob::ODEProblem)
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

    # Update the additional values stored in the tracker
    update_peramp!(tracker, p, (period, amplitude))

    return -fitness
end



"""Custom data structure to store the period and amplitude of each individual"""
struct PeriodAmplitudes
    peramps::Dict{Vector{Float64}, Tuple{Float64, Float64}}

    PeriodAmplitudes() = new(Dict{Vector{Float64}, Tuple{Float64, Float64}}())
end    


"""Helper function to update the period and amplitude of an individual, using in place in CostFunction"""
function update_peramp!(tracker::PeriodAmplitudes, parameters::Vector{Float64}, values::Tuple{Float64, Float64})
    tracker.peramps[parameters] = values
end
