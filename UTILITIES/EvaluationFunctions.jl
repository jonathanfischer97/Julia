
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
    #window = max(1,cld(arrLen,window_ratio)) #* window size is 1% of array length, or 1 if array length is less than 100
    sum_std = sum(std(@view fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs) #* sum rolling window of standard deviations

    # return sum_std / length(fft_peakindxs) #* divide by number of peaks to get average std
end 

"""Return normalized FFT of solution vector. Modifies the solution vector in place"""
function getFrequencies(timeseries) 
    abs.(rfft(timeseries)) ./ cld(length(timeseries), 2)
    # return res ./ cld(length(timeseries), 2) #* normalize by length of timeseries
end


"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution)

    Amem_sol = map(sum, sol.u)

    indx_max, vals_max = findextrema(Amem_sol; height = 1e-2, distance = 5)
    indx_min, vals_min = findextrema(Amem_sol; height = 0.0, distance = 5, find_maxima=false)
    return getPerAmp(sol.t, indx_max, vals_max, indx_min, vals_min)
end

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(solt, indx_max::Vector{Int}, vals_max::Vector{Float64}, indx_min::Vector{Int}, vals_min::Vector{Float64})

    #* Calculate amplitudes and periods
    pers = (solt[indx_max[i+1]] - solt[indx_max[i]] for i in 1:(length(indx_max)-1))
    amps = (vals_max[i] - vals_min[i] for i in 1:min(length(indx_max), length(indx_min)))

    return mean(pers), mean(amps) .|> abs #TODO fix this, why is amps empty sometimes
end



"""Normalize a time series to have mean 0 and amplitude 1 before FFT"""
function normalize_time_series!(ts::Vector{Float64})
    mu = mean(ts)
    amplitude = maximum(ts) - minimum(ts)
    ts .= (ts .- mu) ./ amplitude
end


"""Cost function to be plugged into eval_fitness wrapper"""
function CostFunction(sol::ODESolution)
    Amem_sol = map(sum, sol.u)
    CostFunction(Amem_sol, sol.t)
end


function CostFunction(solu::Vector{Float64}, solt::Vector{Float64})

    tstart = cld(length(solt),10) 

    #* Check if last half of the solution array is steady state
    testwindow = @view solu[end-tstart:end]
    if std(testwindow; mean=mean(testwindow)) < 0.01  
        return [0.0, 0.0, 0.0]
    end 

    #* Trim first 10% of the solution array to avoid initial spikes
    # solu = solu[tstart:end] 

    indx_max, vals_max = findextrema(solu; height = 1e-2, distance = 5)
    indx_min, vals_min = findextrema(solu; height = 0.0, distance = 5, find_maxima=false)
    # indx_min, vals_min = find_minima_between_maxima(solu, indx_max)

    if length(indx_max) < 2 || length(indx_min) < 2 #* if there is no signal in the time domain,
        return [0.0, 0.0, 0.0]
    end
    
    #* Get the rfft of the solution
    fftData = getFrequencies(solu) |> normalize_time_series!


    fft_peakindexes, fft_peakvals = findextrema(fftData; height = 1e-2, distance = 2) #* get the indexes of the peaks in the fft
    # @info length(fft_peakindexes)
    if length(fft_peakindexes) < 2 #* if there is no signal in the frequency domain, return 0.0s
        return [0.0, 0.0, 0.0]
    else
        standard_deviation = getSTD(fft_peakindexes, fftData) #* get the summed standard deviation of the peaks in frequency domain
        sum_diff = getDif(fft_peakvals) #* get the summed difference between the first and last peaks in frequency domain
    
        #* Compute the period and amplitude
        # tview = @view solt[tstart:end]
        period, amplitude = getPerAmp(solt, indx_max, vals_max, indx_min, vals_min)
    
        return [standard_deviation + sum_diff + log(10,period), period, amplitude]
    end
end


#! EVALUATION FUNCTIONS ## 
"""Evaluate the fitness of an individual with new parameters"""
function eval_param_fitness(params::Vector{Float64},  prob::OT; idx::Vector{Int} = [6, 9, 10, 11, 12, 15, 16]) where OT <: ODEProblem
    #* remake with new parameters
    new_prob = remake(prob, p=params)
    return solve_for_fitness_peramp(new_prob, idx)
end

"""Evaluate the fitness of an individual with new initial conditions"""
function eval_ic_fitness(initial_conditions::Vector{Float64}, prob::OT; idx::Vector{Int} = [6, 9, 10, 11, 12, 15, 16]) where OT <: ODEProblem
    #* remake with new initial conditions
    new_prob = remake(prob, u0=initial_conditions)
    return solve_for_fitness_peramp(new_prob, idx)
end

"""Evaluate the fitness of an individual with new initial conditions and new parameters"""
function eval_all_fitness(inputs::Vector{Float64}, prob::OT; idx::Vector{Int} = [6, 9, 10, 11, 12, 15, 16]) where OT <: ODEProblem
    #* remake with new initial conditions
    # new_prob = remake(prob, p = inputs[1:13], u0= [inputs[14:end]; zeros(12)])
    newp = @view inputs[1:13]
    newu = @view inputs[14:end]
    first4u = @view inputs[14:17]
    tend = calculate_tspan(newp, first4u)
    new_prob = remake(prob; p = newp, u0= newu, tspan=(0.0, tend))
    return solve_for_fitness_peramp(new_prob, idx)
end



"""Utility function to solve the ODE and return the fitness and period/amplitude"""
function solve_for_fitness_peramp(prob::OT, idx) where OT <: ODEProblem
    tstart = prob.tspan[2] / 10

    savepoints = tstart:0.1:prob.tspan[2]
    sol = solve(prob, AutoTsit5(Rodas5P()), saveat=savepoints, save_idxs=idx, verbose=false, maxiters=1e6)

    if sol.retcode == ReturnCode.Success
        return CostFunction(sol)
    else
        return [0.0, 0.0, 0.0]
    end
end

"""Calculate tspan based on the slowest reaction rate. 
Simply the reciprocal of the minimum first order rate constants, or the reciprocal of the minimum second order rate constants multiplied by the minimum concentration of the reactants
"""
function calculate_tspan(params, initial_conditions; max_t = 1e5)
    #* Get the minimum rate constant
    min_k, min_k_idx = findmin(params)

    if min_k_idx in (1,4,6,8,10) #* If the minimum rate constant is a second order rate constant, multiply by the minimum concentration of the reactants
        #* Get the minimum concentration of the reactants
        min_conc = minimum(initial_conditions)

        #* Calculate the tspan
        return min(max(10.0, (min_k * min_conc)^-1), max_t)
    else #* If the minimum rate constant is a first order rate constant, simply take the reciprocal
        return min(max(10.0, min_k^-1), max_t)
    end
end

#>MODULE END

#< CUSTOM PEAK FINDER

function findextrema(x::Vector{T}; 
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

# function findextrema(x::AbstractVector{T}; 
#                     height::Union{Nothing,<:Real,NTuple{2,<:Real}}=nothing,
#                     distance::Union{Nothing,Int}=nothing) where {T<:Real}
#     midpts = Vector{Int}(undef, 0)
#     i = 2
#     imax = length(x)

#     while i < imax
#         if x[i-1] < x[i]
#             iahead = i + 1
#             while (iahead < imax) && (x[iahead] == x[i])
#                 iahead += 1
#             end
#             if x[iahead] < x[i]
#                 push!(midpts, (i + iahead - 1) ÷ 2)
#                 i = iahead
#             end
#         end
#         i += 1
#     end 

#     #* Filter by height if needed
#     if !isnothing(height)
#         hmin, hmax = height isa Number ? (height, nothing) : height
#         keepheight = (hmin === nothing || x[midpts] .>= hmin) .& (hmax === nothing || x[midpts] .<= hmax)
#         midpts = midpts[keepheight]
#     end

#     #* Filter by distance if needed
#     if !isnothing(distance)
#         priority = x[midpts]
#         keep = selectbypeakdistance(midpts, priority, distance)
#         midpts = midpts[keep]
#     end

#     extrema_indices = midpts
#     extrema_heights = x[extrema_indices]

#     extrema_indices, extrema_heights
# end

# function find_minima_between_maxima(x::AbstractVector{T}, maxima_indices::Vector{Int}) where {T<:Real}
#     minima_indices = Vector{Int}(undef, length(maxima_indices) - 1)
#     minima_values = Vector{T}(undef, length(maxima_indices) - 1)

#     for i in 1:(length(maxima_indices) - 1)
#         range = maxima_indices[i]:maxima_indices[i+1]
#         minima_indices[i] = argmin(x[range]) + range.start - 1
#         minima_values[i] = x[minima_indices[i]]
#     end

#     return minima_indices, minima_values
# end

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














