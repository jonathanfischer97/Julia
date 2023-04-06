using Plots 
using Catalyst
using DifferentialEquations
# using Peaks
using Statistics
# using BenchmarkTools
# using Logging

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

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
p = [0.125, 411.43070022437774, 42.44753785470635, 57.56861328754449, 1.0, 0.75, 0.001, 0.126, 
    102.45983604574394, 123.39827909372974, 25.805378439197266, 470.69414436040074, 3718.0197650684563]

#initial condition list
u0 = [.01, 3.0, .01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.6, 0.3, 0.2]
L, Lp, LpA, LpAP, LK, LpAK, LpAKL, LpP, LpAPLp, A, K, P = u0

umap = [:L => 0.0, :Lp => 3.0, :LpA => 0.0, :LpAP => 0.0, :LK => 0.0, :LpAK => 0.0, :LpAKL => 0.0, :LpP => 0.0, :LpAPLp => 0.0, :A => 0.9, :K => 0.2, :P => 0.3, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]

#conserved quantities
# ktot=K+LK+LpAK+LpAKL;
# ptot=P+LpP+LpAP+LpAPLp;
# atot=A+LpA+LpAK+LpAKL+LpAP+LpAPLp;
# ltot=L+LK+LpAKL;
# lptot=Lp+LpA+LpAK+LpAKL+LpP+LpAP+2*LpAPLp;
# liptot=ltot+lptot;
# tots = [ktot, ptot, atot, liptot]

#timespan for integration
tspan = (0., 50.)

#solve the reduced ODEs
prob = ODEProblem(fullrn, umap, tspan, p)

sol = solve(prob) #solve adaptively


#plot the results
default(lw = 2, size = (1200, 800))
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

function oscillatory_strength(sol; min_freq=0.02, max_freq=0.5)
    # Extract the time points and solution values
    time_points = sol.t
    solution = sol[1,:]

    # Calculate the sampling rate
    sampling_rate = 1 / mean(diff(time_points))

    # Compute the FFT of the signal
    fft_result = rfft(solution)

    # Calculate the power spectrum
    power_spectrum = abs.(fft_result).^2

    # Convert frequency range to indices
    min_freq_index = max(2, round(Int, min_freq * length(solution) / sampling_rate) + 1)
    max_freq_index = min(div(length(power_spectrum), 2), round(Int, max_freq * length(solution) / sampling_rate) + 1)

    # Find the dominant frequency and its power within the specified range
    max_power_index = argmax(power_spectrum[min_freq_index:max_freq_index]) + min_freq_index - 1
    max_power = power_spectrum[max_power_index]

    # Calculate the total power (excluding the DC component)
    total_power = sum(power_spectrum) - power_spectrum[1]

    # Return the ratio of the dominant frequency's power to the total power
    return max_power / total_power
end



function eval_fitness2(p::Vector{Float64},  prob::ODEProblem)
    Y = nothing
    try 
        Y = solve(remake(prob, p=p, saveat=0.01, maxiters=10000, verbose=false))
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
    return -oscillatory_strength(Y[1,:])
end




function make_fitness_function(func::Function, prob::ODEProblem) # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return func(p, prob)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(eval_fitness2, prob) # Create a fitness function that includes your ODE problem as a constant



## parameter constraint ranges ##
ka_min, ka_max = 0.001, 1.
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
    "y" => Dict("min" => 100., "max" => 2000.)
);

random_p = [rand(param_values[p]["min"]:0.01:param_values[p]["max"]) for p in keys(param_values)]
eval_fitness2(random_p, tots, prob)
# eval_fitness(random_p, tots, prob2, 1)
testprob = ODEProblem(reduced_oscillator_odes!, u0[1:2:3], tspan, vcat(random_p, tots))
testsol = solve(testprob, p=vcat(p, tots))
findmaxima(testsol[1,:], 10)
plot(testsol)

#Optimization parameters
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=10, parallelization=:thread, abstol=1e-2, reltol=1e-2)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 1000, selection = tournament(100),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)


# function callback(trace)

#Optimization
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(prob, p=newp))


#plot the results
plot(newsol, xlabel="Time (s)", ylabel="Concentration (mM)", title="Optimized Model")

eval_fitness2(newp, prob)

eval_fitness(newp, tots, prob, 1)
eval_fitness(newp, tots, prob2, 1)


any(isless.(newsol[1,:], 0.0))


newsolfft = abs.(rfft(newsol[1,:]))
plot(newsolfft[4:end])



oscillatory_strength(newsol[1,:])