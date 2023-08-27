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
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData::Vector{Float64}; window::Int =1) #get average standard deviation of fft peak indexes
    arrLen = length(fft_arrayData)
    # window = max(1,cld(arrLen,window_ratio)) #* window size is 1% of array length, or 1 if array length is less than 100
    # window = 5
    sum_std = sum(std(fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations
    return sum_std
    # return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end 

"""Return normalized FFT of solution vector. Modifies the solution vector in place"""
function getFrequencies(timeseries) 
    res = abs.(rfft(timeseries))
    return res ./ cld(length(timeseries), 2) #* normalize by length of timeseries
end

function flip_about_mean(vec::Vector{Float64})
    mean_value = mean(vec)
    return [2 * mean_value - x for x in vec]
end

function flip_about_mean!(vec::Vector{Float64})
    mean_value = mean(vec)
    @inbounds for i in eachindex(vec)
        vec[i] = 2 * mean_value - vec[i]
    end
    vec
end



"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, idx::Int = 1)

    indx_max, vals_max = findextrema(sol[idx,:]; height = 1e-2, distance = 10)
    return getPerAmp(sol, indx_max, vals_max, idx)
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, indx_max::Vector{Int}, vals_max::Vector{Float64}, idx::Int = 1)
    #* Find peaks of the minima too 
    indx_min, vals_min = findextrema(sol[idx,:]; height = 1e-2, distance = 10, find_maxima=false)

    #* Calculate amplitudes and periods
    # vals_max = maxprops["peak_heights"]
    # vals_min = minprops["peak_heights"]

    #* Calculate amplitudes and periods

    pers = (sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1))
    amps = ((vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min)))


    return maximum(pers), mean(amps) .|> abs #TODO fix this, why is amps empty sometimes
end


"""Normalize a time series to have mean 0 and amplitude 1 before FFT"""
function normalize_time_series!(ts::Vector{Float64})::Vector{Float64}
    mu = mean(ts)
    amplitude = maximum(ts) - minimum(ts)
    ts .= (ts .- mu) ./ amplitude
end


"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution; idx::Int = 1)::Vector{Float64}

    #* Check if last half of the solution array is steady state
    lasthalfsol = sol[cld(length(sol.t),2):end]
    if std(lasthalfsol[idx,:]) < 0.05 
        return [0.0, 0.0, 0.0]
    end 

    #* Trim first 10% of the solution array to avoid initial spikes
    tstart = cld(length(sol.t),10) 
    trimsol = sol[tstart:end] 

    indx_max, vals_max = findextrema(trimsol[idx,:]; height = 1e-2, distance = 10)
    indx_min, vals_min = findextrema(trimsol[idx,:]; height = 1e-2, distance = 10, find_maxima=false)

    if length(indx_max) < 2 || length(indx_min) < 2 #* if there is no signal in the frequency domain,
        return [0.0, 0.0, 0.0]
    end
    
    #* Get the rfft of the solution
    fftData = getFrequencies(trimsol[idx,:])

    #* Normalize the solution array. WARNING: solarray is modified after this line
    normalize_time_series!(fftData)

    fft_peakindexes, fft_peakvals = findextrema(fftData; height = 1e-3, distance = 2) #* get the indexes of the peaks in the fft
    # @info length(fft_peakindexes)
    if length(fft_peakindexes) < 2 #* if there is no signal in the frequency domain, return 0.0s
        return [0.0, 0.0, 0.0]
    else
        # fft_peakvals = peakprops["peak_heights"]

        standard_deviation = getSTD(fft_peakindexes, fftData; window = 5) #* get the summed standard deviation of the peaks in frequency domain
        sum_diff = getDif(fft_peakvals) #* get the summed difference between the first and last peaks in frequency domain
    
        #* Compute the period and amplitude
        period, amplitude = getPerAmp(sol, indx_max, vals_max, idx)
    
        return [-standard_deviation - sum_diff - log(10,period), period, amplitude]
    end
end


#! EVALUATION FUNCTIONS ## 
"""Evaluate the fitness of an individual with new parameters"""
function eval_param_fitness(params::Vector{Float64},  prob::ODEProblem; idx::Int = 4)
    #* remake with new parameters
    new_prob = remake(prob, p=params)
    return solve_for_fitness_peramp(new_prob, idx)
end

"""Evaluate the fitness of an individual with new initial conditions"""
function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::ODEProblem; idx::Int = 4)
    #* remake with new initial conditions
    new_prob = remake(prob, u0=[initial_conditions; zeros(length(prob.u0)-length(initial_conditions))])
    return solve_for_fitness_peramp(new_prob, idx)
end

"""Evaluate the fitness of an individual with new initial conditions and new parameters"""
function eval_all_fitness(params::Vector{Float64}, initial_conditions::Vector{Float64}, prob::ODEProblem; idx::Int = 4)
    #* remake with new initial conditions
    new_prob = remake(prob, p = params, u0=[initial_conditions; zeros(length(prob.u0)-length(initial_conditions))])
    return solve_for_fitness_peramp(new_prob, idx)
end

"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function solve_for_fitness_peramp(prob::ODEProblem, idx::Int = 4)

    sol = solve(prob, Rosenbrock23(), saveat=0.1, save_idxs=idx, verbose=false)
    # return CostFunction(sol)
    
    if sol.retcode == ReturnCode.Success
        return CostFunction(sol)
    else
        return [0.0, 0.0, 0.0]
    end
end

#>MODULE END

#< CUSTOM PEAK FINDER

function findextrema(x::AbstractVector{T}; 
                    height::Union{Nothing,<:Real,NTuple{2,<:Real}}=nothing,
                    distance::Union{Nothing,Int}=nothing,
                    find_maxima::Bool=true) where {T<:Real}
    midpts = Vector{Int}(undef, 0)
    i = 2
    imax = length(x)

    while i < imax
        if (find_maxima && x[i-1] < x[i]) || (!find_maxima && x[i-1] > x[i])
            iahead = i + 1
            while (iahead < imax) && (x[iahead] == x[i])
                iahead += 1
            end

            if (find_maxima && x[iahead] < x[i]) || (!find_maxima && x[iahead] > x[i])
                push!(midpts, (i + iahead - 1) ÷ 2)
                i = iahead
            end
        end
        i += 1
    end 

    #* Filter by height if needed
    if !isnothing(height)
        hmin, hmax = height isa Number ? (height, nothing) : height
        keepheight = (hmin === nothing || x[midpts] .>= hmin) .& (hmax === nothing || x[midpts] .<= hmax)
        midpts = midpts[keepheight]
    end

    #* Filter by distance if needed
    if !isnothing(distance)
        priority = find_maxima ? x[midpts] : -x[midpts] # Use negative values for minima
        keep = selectbypeakdistance(midpts, priority, distance)
        midpts = midpts[keep]
    end

    extrema_indices = midpts
    extrema_heights = x[extrema_indices]

    extrema_indices, extrema_heights
end

function selectbypeakdistance(pkindices, priority, distance)
    npkindices = length(pkindices)
    keep = trues(npkindices)

    prioritytoposition = sortperm(priority, rev=true)
    for i ∈ npkindices:-1:1
        j = prioritytoposition[i]
        iszero(keep[j]) && continue

        k = j-1
        while (1 <= k) && ((pkindices[j]-pkindices[k]) < distance)
            keep[k] = false
            k -= 1
        end

        k = j+1
        while (k <= npkindices) && ((pkindices[k]-pkindices[j]) < distance)
            keep[k] = false
            k += 1
        end
    end
    keep
end

# """
# Struct to hold the properties of the peaks found by the peak finder
# """
# struct PeakProperties
#     peak_heights::Union{Nothing, Vector{Float64}}
#     prominences::Union{Nothing, Vector{Float64}}
#     leftbases::Union{Nothing, Vector{Int}}
#     rightbases::Union{Nothing, Vector{Int}}
#     widths::Union{Nothing, Vector{Float64}}
#     widthheights::Union{Nothing, Vector{Float64}}
#     leftips::Union{Nothing, Vector{Float64}}
#     rightips::Union{Nothing, Vector{Float64}}
# end


# function filterproperties!(properties::PeakProperties, keep::BitVector)
#     properties.peak_heights = properties.peak_heights[keep]
#     properties.prominences = properties.prominences[keep]
#     properties.leftbases = properties.leftbases[keep]
#     properties.rightbases = properties.rightbases[keep]
#     properties.widths = properties.widths[keep]
#     properties.widthheights = properties.widthheights[keep]
#     properties.leftips = properties.leftips[keep]
#     properties.rightips = properties.rightips[keep]
# end

# function findpeaks1d(x::AbstractVector{T};
#                     height::Union{Nothing,<:Real,NTuple{2,<:Real}}=nothing,
#                     distance::Union{Nothing,I}=nothing,
#                     prominence::Union{Nothing,Real,NTuple{2,Real}}=nothing,
#                     width::Union{Nothing,Real,NTuple{2,Real}}=nothing,
#                     wlen::Union{Nothing,I}=nothing,
#                     relheight::Real=0.5,
#                     calc_peak_heights::Bool=false,
#                     calc_prominences::Bool=false,
#                     calc_widths::Bool=false) where {T<:Real,I<:Integer}

#     pkindices, leftedges, rightedges = localmaxima1d(x)

#     # Initialize variables for optional calculations
#     peak_heights = nothing
#     prominences = nothing
#     leftbases = nothing
#     rightbases = nothing
#     widths = nothing
#     widthheights = nothing
#     leftips = nothing
#     rightips = nothing

#     if calc_peak_heights && !isnothing(height)
#         pkheights = x[pkindices]
#         hmin, hmax = height isa Number ? (height, nothing) : height
#         keepheight = selectbyproperty(pkheights, hmin, hmax)
#         pkindices = pkindices[keepheight]
#         peak_heights = pkheights[keepheight]
#     end

#     if !isnothing(distance)
#         keepdist = selectbypeakdistance(pkindices, x[pkindices], distance)
#         pkindices = pkindices[keepdist]
#     end

#     if calc_prominences && (!isnothing(prominence) || !isnothing(width))
#         prominences, leftbases, rightbases = peakprominences1d(x, pkindices, wlen)
#     end

#     if !isnothing(prominence)
#         pmin, pmax = prominence isa Number ? (prominence, nothing) : prominence
#         keepprom = selectbyproperty(prominences, pmin, pmax)
#         pkindices = pkindices[keepprom]
#     end

#     if calc_widths && !isnothing(width)
#         widths, widthheights, leftips, rightips = peakwidths1d(x, pkindices, relheight, prominences, leftbases, rightbases)
#         wmin, wmax = width isa Number ? (width, nothing) : width
#         keepwidth = selectbyproperty(widths, wmin, wmax)
#         pkindices = pkindices[keepwidth]
#     end

#     # Construct the properties struct with the calculated values
#     properties = PeakProperties(peak_heights, prominences, leftbases, rightbases, widths, widthheights, leftips, rightips)

#     pkindices, properties
# end














