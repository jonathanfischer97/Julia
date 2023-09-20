using Plots 
using Catalyst
using DifferentialEquations
using Statistics
using FFTW
using Peaks


"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*y,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 

    #previously hidden NERDSS reactions
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*y,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*y,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*y,kb3), LpA + LK <--> LpAKL
    (ka4*y,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end  

ogrn = @reaction_network ogrn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*y,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 
end 

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
p = [0.125, 411.43070022437774, 42.44753785470635, 57.56861328754449, 1.0, 0.75, 0.001, 0.126, 
    102.45983604574394, 123.39827909372974, 25.805378439197266, 470.69414436040074, 3718.0197650684563]

pmap = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848,
    :ka2 => 1.0, :kb2 => 0.9504699043910143,
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123,
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.]

#initial condition list
umap = [:L => 0.0, :Lp => 3.0, :LpA => 0.0, :LpAP => 0.0, :LK => 0.0, :LpAK => 0.0, :LpAKL => 0.0, :LpP => 0.0, :LpAPLp => 0.0, :A => 0.9, :K => 0.2, :P => 0.3, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]


#timespan for integration
tspan = (0., 100.)

#solve ODEs for the full model
fullprob = ODEProblem(fullrn, umap, tspan, pmap)
#solve ODEs for the original model
ogprob = ODEProblem(ogrn, umap[1:end-4], tspan, pmap)

fullsol = solve(fullprob, saveat = 0.01) 
ogsol = solve(ogprob, saveat = 0.01) 


#plot the results
default(lw = 2, size = (1200, 800))
plot(fullsol)
plot!(ogsol, ls = :dash, legend = false)




## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get summed difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = sum(diff(arrayData[indexes])) + arrayData[indexes[end]] #add the last element
    return sum_diff
end

function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
    window = max(1, round(Int, window_ratio * length(arrayData)))
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in peakindxs)
    return sum_std / length(peakindxs)
end

function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    return res ./ cld(length(y), 2) #normalize amplitudes
end

function eval_fitness(p::Vector{Float64},  prob::ODEProblem)
    Y = solve(remake(prob, p=p, saveat=0.01))
    #get the fft of the solution
    fftData = getFrequencies(Y[1,:])
    indexes = findmaxima(fftData,10)[1] #get the indexes of the peaks in the fft
    timepeaks = length(findmaxima(fftData,1)[1]) #get the times of the peaks in the fft
    if isempty(indexes) || timepeaks < 2 #if there are no peaks, return 0
        return 0.0
    end
    std = getSTD(indexes, fftData, 0.1) #get the standard deviation of the peaks
    diff = getDif(indexes, fftData) #get the difference between the peaks
    return -std - diff
end

"""Eval fitness function for testing"""
function eval_fitness(t, sol::Vector, target::Float64, bias::Float64)
    #get the fft of the solution
    fftData = getFrequencies(sol)
    indexes = findmaxima(fftData,1)[1] #get the indexes of the peaks in the fft
    #timepeaks = length(findmaxima(fftData,1)[1]) #get the times of the peaks in the fft
    if isempty(indexes) #|| timepeaks < 2 #if there are no peaks, return 0
        return 0.0
    end
    std = getSTD(indexes, fftData, 0.1) #get the standard deviation of the peaks
    @info "std: $std"
    diff = getDif(indexes, fftData) #get the difference between the peaks
    @info "diff: $diff"
    return std + diff
end



"""Test cost function"""
function oscillatory_strength(t, sol::Vector; min_freq=0.02, max_freq=0.5)
    # Extract the time points and solution values
    time_points = t
    solution = sol

    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the power spectrum
    power_spectrum = abs.(fft_result).^2

    # Convert frequency range to indices
    min_freq_index = max(2, round(Int, min_freq * length(signal) / sampling_rate) + 1)
    max_freq_index = min(div(length(power_spectrum), 2), round(Int, max_freq * length(signal) / sampling_rate) + 1)

    # Find the dominant frequency and its power within the specified range
    max_power_index = argmax(power_spectrum[min_freq_index:max_freq_index]) + min_freq_index - 1
    max_power = power_spectrum[max_power_index]

    # Calculate the total power (excluding the DC component)
    total_power = sum(power_spectrum) - power_spectrum[1]

    # Return the ratio of the dominant frequency's power to the total power
    return max_power / total_power
end

"""Amplitude spectrum analysis, biased to target frequency, exponential decay"""
function oscillatory_strength_biased(t, sol::Vector, target_frequency::Float64, bias_strength::Float64=0.1)
    solution = sol
    time_points = t
    
    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the amplitude spectrum
    amplitude_spectrum = abs.(fft_result)

    # Calculate frequency resolution
    frequency_resolution = sampling_rate / length(solution)

    # Calculate the index of the target frequency
    target_index = round(Int, target_frequency / frequency_resolution) + 1

    # Compute a weight function that biases the cost towards the target frequency
    weight_function = [exp(-bias_strength * abs(i - target_index)) for i in 1:length(fft_result)]

    # Calculate the weighted amplitude
    weighted_amplitude_spectrum = weight_function .* amplitude_spectrum

    # Find the dominant frequency and its weighted amplitude
    max_weighted_amplitude_index = argmax(weighted_amplitude_spectrum[2:end]) + 1
    max_weighted_amplitude = weighted_amplitude_spectrum[max_weighted_amplitude_index]

    # Calculate the total weighted amplitude (excluding the DC component)
    total_weighted_amplitude = sum(weighted_amplitude_spectrum) - weighted_amplitude_spectrum[1]

    # Return the ratio of the dominant frequency's weighted amplitude to the total weighted amplitude
    return max_weighted_amplitude / total_weighted_amplitude
end

"""Unbiased amplitude spectrum analysis, exponential decay"""
function oscillatory_strength_unbiased(t, sol::Vector, target_frequency::Float64, bias_strength::Float64=0.1)
    solution = sol
    time_points = t
    
    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the amplitude spectrum
    amplitude_spectrum = abs.(fft_result)

    # Calculate frequency resolution
    frequency_resolution = sampling_rate / length(solution)

    # Find the dominant frequency and its amplitude
    max_amplitude_index = argmax(amplitude_spectrum[2:end]) + 1
    max_amplitude = amplitude_spectrum[max_amplitude_index]

    # Calculate the total amplitude (excluding the DC component)
    total_amplitude = sum(amplitude_spectrum) - amplitude_spectrum[1]

    # Return the ratio of the dominant frequency's amplitude to the total amplitude
    return max_amplitude / total_amplitude
end

"""Calculates the target power over the weighted total power"""
function oscillatory_strength_sigmoid(t, sol::Vector, target_freq::Float64, bias_sigma::Float64=0.1)
    solution = sol
    time_points = t
    
    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the power spectrum
    power_spectrum = abs.(fft_result).^2

    # Convert the target frequency to an index
    target_freq_index = max(2, round(Int, target_freq * length(signal) / sampling_rate) + 1)

    # Create a double sigmoid weight function centered at the target frequency index
    freq_indices = 2:length(power_spectrum)
    double_sigmoid_weights = [1 / (1 + exp(-bias_sigma * (i - target_freq_index))) * 1 / (1 + exp(bias_sigma * (i - target_freq_index))) for i in freq_indices]

    # Calculate the weighted total power (excluding the DC component)
    weighted_total_power = sum(double_sigmoid_weights .* power_spectrum[freq_indices])

    # Make sure the weighted total power is not zero to avoid division by zero
    if weighted_total_power == 0.0
        return 0.0
    end

    # Calculate the power at the target frequency index
    target_power = power_spectrum[target_freq_index]

    # Return the ratio of the target frequency's power to the weighted total power
    return target_power / weighted_total_power
end







# Test with a sine wave (oscillatory signal)
t = collect(range(0, stop=10π, length=10001))



#Plot score vs frequency for biased cost function
function get_sincosts(target, bias, costfunc)
    freqs = collect(range(0.01, stop=1., length=100))
    scores = zeros(100)
    for i in 1:100
        sine_wave = sin.(2π * freqs[i] * t)
        scores[i] = costfunc(t, sine_wave, target, bias)
    end
    return freqs, scores
end

function get_randcosts(target, bias, costfunc)
    freqs = collect(range(0.01, stop=1., length=100))
    scores = zeros(100)
    for i in 1:100
        random_signal = rand(10001)
        scores[i] = costfunc(t, random_signal, target, bias)
    end
    return freqs, scores
end

function get_flatcosts(target, bias, costfunc)
    freqs = collect(range(0.01, stop=1., length=100))
    scores = zeros(100)
    for i in 1:100
        flat_signal = ones(10001)
        scores[i] = costfunc(t, flat_signal, target, bias)
    end
    return freqs, scores
end

function plotcosts(costfunc, target = 0.7)
    minbias = get_sincosts(target, 0.1, costfunc)
    midbias = get_sincosts(target, 0.5, costfunc)
    maxbias = get_sincosts(target, 1.9, costfunc)
    randminbias = get_randcosts(target, 0.1, costfunc)
    flatminbias = get_flatcosts(target, 0.1, costfunc)

    costplot = plot(minbias, title="Target: $target", label="Bias: 0.1", xlabel="Frequency", ylabel="Score")
    plot!(costplot, midbias, label="Bias: 0.5")
    plot!(costplot, maxbias, label="Bias: 0.9")
    plot!(costplot, randminbias, label="Random signal, bias: 0.1")
    plot!(costplot, flatminbias, label="Flat signal, bias: 0.1")
    # scatter!(costplot, [target], [1], label="Target frequency")
    display(costplot)
end

plotcosts(oscillatory_strength_unbiased, 0.5)

#double barplot of score of fullsol and ogsol for different cost functions
bar([eval_fitness(t, fullsol, 0.5, 0.1), eval_fitness(t, ogsol, 0.5, 0.1)], label=["Full" "OG"], title="Score of fullsol and ogsol for different cost functions", ylabel="Score", xticks=(1:2, ["Oscillatory" "Random"]), yticks=0:0.1:1, ylims=(0,1), legend=:topleft)


pmapcopy = deepcopy(pmap)
pcopy = [x[2] for x in pmapcopy]
fullcosts = []
ogcosts = []
prange = 100.:100.:5000
for y in prange
    # update parm map with new y
    pmapcopy[end] = :y => y
    # println("pmapcopy: $pmapcopy")
    fullprob = ODEProblem(fullrn, umap, tspan, pmapcopy)
    ogprob = ODEProblem(ogrn, umap[1:end-4], tspan, pmapcopy)

    
    # solve both full and original model with new y
    fullsol = solve(fullprob, saveat=0.01)
    ogsol = solve(ogprob, saveat=0.01)

    # evaluate fitness of both solutions
    fullcost = oscillatory_strength_unbiased(t, fullsol[1,:], 0.5, 0.1)
    ogcost = oscillatory_strength_unbiased(t, ogsol[1,:], 0.5, 0.1)

    # push to array
    push!(fullcosts, fullcost)
    push!(ogcosts, ogcost)
end
# plot
p1 = plot(prange, fullcosts, label="Full", title="Score of fullsol and ogsol for different y", ylabel="Score", xlabel="y", legend=:topleft)
plot!(p1, prange, ogcosts, label="OG")
display(p1)


pmapcopy[end] = :y => 5000.
ogprob = ODEProblem(ogrn, umap[1:end-4], tspan, pmapcopy)
testsol = solve(ogprob, saveat=0.01)
plot(testsol, title="Full model", xlabel="Time", ylabel="x, y", legend=:topleft)
eval_fitness(t, testsol[1,:], 0.5, 0.1)