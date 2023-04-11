using Plots 
using Catalyst
using DifferentialEquations
using Statistics
default(lw = 2, size = (1200, 800))

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
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

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
pmap = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848,
    :ka2 => 1.0, :kb2 => 0.9504699043910143,
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123,
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 2000.]
p = [x[2] for x in pmap]
    
#initial condition list
umap = [:L => 0.0, :Lp => 3.0, :K => 0.2, :P => 0.3, :A => 0.9, :LpA => 0.0, :LK => 0.0, :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in umap]

#timespan for integration
tspan = (0., 100.)

#solve the reduced ODEs
prob = ODEProblem(fullrn, umap, tspan, p)

sol = solve(prob) #solve adaptively

#plot the results
plot(sol)





## Genetic algorithm to find the best parameters for the reduced model ## 
using Evolutionary
using FFTW

## COST FUNCTIONS and dependencies 
# ChatGPT defined replacement function
function oscillatory_strength(solution::Vector)
    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the power spectrum
    power_spectrum = abs.(fft_result).^2

    # Find the dominant frequency and its power
    max_power_index = argmax(power_spectrum[2:div(end, 2)]) + 1
    max_power = power_spectrum[max_power_index]

    # Calculate the total power (excluding the DC component)
    total_power = sum(power_spectrum) - power_spectrum[1]

    # Return the ratio of the dominant frequency's power to the total power
    return max_power / total_power
end


"""Unbiased amplitude spectrum analysis, exponential decay"""
function oscillatory_strength_unbiased(sol::ODESolution)
    solution = sol[1,:]
    time_points = sol.t
    
    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Remove the mean of the signal to eliminate the DC component
    signal = solution .- mean(solution)

    # Compute the FFT of the signal
    fft_result = rfft(signal)

    # Calculate the amplitude spectrum
    amplitude_spectrum = abs.(fft_result)

    # Calculate frequency resolution
    # frequency_resolution = sampling_rate / length(solution)

    # Find the dominant frequency and its amplitude
    max_amplitude_index = argmax(amplitude_spectrum[2:end]) + 1
    max_amplitude = amplitude_spectrum[max_amplitude_index]

    # Calculate the total amplitude (excluding the DC component)
    total_amplitude = sum(amplitude_spectrum) - amplitude_spectrum[1]

    # Return the ratio of the dominant frequency's amplitude to the total amplitude
    return max_amplitude / total_amplitude
end







"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(p::Vector{Float64},  prob::ODEProblem)
    Y = nothing
    pvals = p[1:13]
    uvals = p[14:end]
    try 
        Y = solve(remake(prob, p=pvals, u0=vcat(uvals, fill(0.0,11))), saveat=0.01, maxiters=10000, verbose=false)
        if Y.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(Y) for x in array) || any(x==1 for array in isless.(Y, 0.0) for x in array)
            return 1.0
        end
    catch e 
        if isa(e, DomainError) #catch domain errors
            return 1.0
        else
            rethrow(e) #rethrow other errors
        end
    end
    return -oscillatory_strength_unbiased(Y)
end


"""Regular cost function"""
function eval_fitness(p::Vector{Float64},  prob::ODEProblem)
    Y = solve(remake(prob, p=p), saveat=0.01, save_idxs=1)
    return -oscillatory_strength_unbiased(Y)
end




function make_fitness_function(func::Function, prob::ODEProblem) # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return func(p, prob)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(eval_fitness_catcherrors, prob) # Create a fitness function that includes your ODE problem as a constant



## parameter constraint ranges ##
ka_min, ka_max = 0.001, 10.
kb_min, kb_max = 0.001, 500.0
kcat_min, kcat_max = 0.001, 500.0


param_values = Dict(
    "ka1" => Dict("min" => ka_min, "max" => ka_max),
    "kb1" => Dict("min" => kb_min, "max" => kb_max),
    "kcat1" => Dict("min" => kcat_min, "max" => kcat_max),
    "ka2" => Dict("min" => ka_min, "max" => ka_max),
    "kb2" => Dict("min" => kb_min, "max" => kb_max),
    "ka3" => Dict("min" => ka_min, "max" => ka_max),
    "kb3" => Dict("min" => kb_min, "max" => kb_max),
    "ka4" => Dict("min" => ka_min, "max" => ka_max),
    "kb4" => Dict("min" => kb_min, "max" => kb_max),
    "ka7" => Dict("min" => ka_min, "max" => ka_max),
    "kb7" => Dict("min" => kb_min, "max" => kb_max),
    "kcat7" => Dict("min" => kcat_min, "max" => kcat_max),
    "y" => Dict("min" => 100., "max" => 5000.),

    "L" => Dict("min" => 0.0, "max" => 3.0),
    "Lp" => Dict("min" => 0.0, "max" => 3.0),
    "K" => Dict("min" => 0.0, "max" => 1.0),
    "P" => Dict("min" => 0.0, "max" => 1.0),
    "A" => Dict("min" => 0.0, "max" => 3.0),
);

#Optimization parameters
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=15, parallelization=:thread, abstol=1e-2, reltol=1e-2)

common_range = 0.5
valrange = fill(common_range, 18)
mthd = GA(populationSize = 5000, selection = tournament(500),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)

#Optimization
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(prob, p=newp))


#plot the results
plot(newsol, xlabel="Time (s)", ylabel="Concentration (mM)", title="Optimized Model")
