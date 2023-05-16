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
    using Unitful: µM, nm, s
    using BenchmarkTools, Profile, ProgressMeter
    # using MultivariateStats, UMAP, TSne, StatsPlots
    using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
end

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase, P = Phosphatase, A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

    #previously excluded reactions, all possible combinations possible in vitro
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

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
#parameter list
    const param_names = ["ka1", "kb1", "kcat1", "ka2", "kb2", "ka3", "kb3", "ka4", "kb4", "ka7", "kb7", "kcat7", "V/A"]

    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :y => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #initial condition list
    usym = [:L => 0.0, :Lp => 3.0, :K => 0.5, :P => 0.3, :A => 2.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #timespan for integration
    tspan = (0., 100.)
    #solve the reduced ODEs
    prob = ODEProblem(fullrn, u0, tspan, p)
    sol = solve(prob, saveat=0.1) #solve adaptively
    #* plot the results save
    plot(sol)
end


#! Helper functions for cost function ## 
begin
    """Get summed difference of peaks in the frequency domain"""
    function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) 
        idxarrLen = length(indexes)
        if idxarrLen < 2
            return 0.0
        end
        sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(idxarrLen-1))
        sum_diff += arrayData[indexes[end]]
        return sum_diff
    end

    """Get summed average standard deviation of peaks in the frequency domain"""
    function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
        arrLen = length(arrayData)
        window = max(1, round(Int, window_ratio * arrLen))
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
            @inbounds amps = [(vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min))]
            @inbounds pers = [sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1)]

            # Calculate means of amplitudes and periods
            amp = mean(amps)
            per = mean(pers)

            return per, amp
        end
    end

    """Cost function to be plugged into eval_fitness wrapper"""
    function CostFunction(Y::ODESolution)
        #get the fft of the solution
        fftData = getFrequencies(Y.u)
        fftindexes = findmaxima(fftData,1)[1] #get the indexes of the peaks in the fft
        timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the fft
        if isempty(fftindexes) || length(timeindexes) < 2 #if there are no peaks, return 0
            return 0.0, 0.0, 0.0
        end
        std = getSTD(fftindexes, fftData, 0.0001) #get the standard deviation of the peaks
        diff = getDif(fftindexes, fftData) #get the difference between the peaks

        # Compute the period and amplitude
        period, amplitude = getPerAmp(Y)

        # Return cost, period, and amplitude as a tuple
        return -std + diff, period, amplitude
    end
end

begin
    """Custom data structure to store the period and amplitude of each individual"""
    struct PeriodAmplitudes
        peramps::Dict{Vector{Float64}, Tuple{Float64, Float64}}
    
        PeriodAmplitudes() = new(Dict{Vector{Float64}, Tuple{Float64, Float64}}())
    end    

    function update_peramp!(tracker::PeriodAmplitudes, parameters::Vector{Float64}, values::Tuple{Float64, Float64})
        tracker.peramps[parameters] = values
    end
end


"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(tracker::PeriodAmplitudes, p::Vector{Float64},  prob::ODEProblem)
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

# function double_eval(eval_function, sol)
#     firstscore = eval_function()
#     secondscore = eval_function(sol[length(sol),:])

#     if firstscore = 0.0 || secondscore = 0.0
#         return 0.0
#     else
#         return firstscore + secondscore
#     end

begin
    ## parameter constraint ranges ##
    ka_min, ka_max = 0.001, 10. #uM^-1s^-1
    kb_min, kb_max = 0.001, 1000.0 #s^-1
    kcat_min, kcat_max = 0.001, 1000.0 #s^-1

    param_values = OrderedDict(
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
        "y" => Dict("min" => 10., "max" => 20000.)
        # "L" => Dict("min" => 0.1, "max" => 3.0),
        # "Lp" => Dict("min" => 0.1, "max" => 3.0),
        # "K" => Dict("min" => 0.1, "max" => 1.0),
        # "P" => Dict("min" => 0.1, "max" => 1.0),
        # "A" => Dict("min" => 0.1, "max" => 3.0),
    );
end

"""Function that generates population of log-uniform sampled random parameter values"""
function generate_population(param_values::OrderedDict, n::Int)
    params = keys(param_values)
    n_params = length(params)
    population = Matrix{Float64}(undef, n_params, n)
    for i in 1:n
        for (ind, param) in enumerate(params)
            min_val = log(param_values[param]["min"])
            max_val = log(param_values[param]["max"])
            min_val < max_val ? population[ind, i] = exp(rand(Uniform(min_val, max_val))) : population[ind, i] = exp(min_val)
        end
    end
    return [collect(population[:, i]) for i in 1:n]
end


#! OVERRIDES FOR Evolutionary.jl ##
"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method::GA, options)
    record["populationmap"] = [(ind=population[i], fit=state.fitpop[i]) for i in eachindex(population)]
end

"""Show override function to prevent printing large arrays"""
function Evolutionary.show(io::IO, t::Evolutionary.OptimizationTraceRecord)
    print(io, lpad("$(t.iteration)",6))
    print(io, "   ")
    print(io, lpad("$(t.value)",14))
    for (key, value) in t.metadata
        if !isa(value, AbstractArray)
            print(io, "\n * $key: $value")
        end
    end
    return
end



#! SINGLE RUN GENETIC ALGORITHM ##
"""Wrapper function to create a fitness function that includes your ODE problem as a constant"""
function make_fitness_function(func::Function, prob::ODEProblem)
    function fitness_function(tracker::PeriodAmplitudes, p::Vector{Float64})
        return func(tracker, p, prob)
    end
    return fitness_function
end

# Create a PeriodAmplitudes instance
tracker = PeriodAmplitudes()

fitness_function = make_fitness_function(eval_fitness_catcherrors, prob) # Create a fitness function that includes your ODE problem as a constant

# Wrap the fitness_function call to work with the Evolutionary.jl optimize function
wrapped_fitness(params) = fitness_function(tracker, params)


#! Optimization block
begin
    population_size = 10000
    pop = generate_population(param_values, population_size)
    constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
    opts = Evolutionary.Options(abstol=1e-2, reltol=1.00, successive_f_tol = 5, iterations=10, store_trace = true, 
            show_trace=true, show_every=1, parallelization=:thread)
    common_range = 0.5; valrange = fill(common_range, length(param_values))
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
            crossover = TPX, crossoverRate = 0.5,
            mutation  = BGA(valrange, 2), mutationRate = 0.7)

    # Progress Bar
    progress = ProgressThresh(opts.abstol, "Converging: ")

    # Optimization
    result = Evolutionary.optimize(wrapped_fitness, constraints, mthd, pop, opts)

    # Get the population map
    fitpops = [gen.metadata["populationmap"] for gen in result.trace]
    fitpops = reduce(vcat, fitpops)
    # Filter out individuals with fitness values less than 0.4
    filtered_fitpops = filter(x -> x.fit < -0.4, fitpops)

    # Get the best solution
    newp = result.minimizer
    newsol = solve(remake(prob; p=newp))

    # Access the additional values after the optimization is finished
    additional_values = tracker.peramps

    #plot the results
    plot(newsol, xlabel="Time (s)", ylabel="Concentration (mM)", title="Optimized Model")
end

fitpops = [gen.metadata["populationmap"] for gen in result.trace]
fitpops = reduce(vcat, fitpops)
# Filter out individuals with fitness values less than 0.4
filtered_fitpops = filter(x -> x.fit < -0.4, fitpops)
fitinds = [x.ind for x in filtered_fitpops]


"""Iterate through the dataset, simulate the system, and calculate the amplitude and period"""
function analyze_population(result::Evolutionary.OptimizationResults)
    fitpops = [gen.metadata["populationmap"] for gen in result.trace]
    fitpops = reduce(vcat, fitpops)
    # Filter out individuals with fitness values less than 0.4
    filtered_fitpops = filter(x -> x.fit < -0.4, fitpops)

    extended_dataset = []

    for individual in filtered_fitpops
        ind = individual.ind
        fitness = individual.fit
        sol = solve(remake(prob; p=ind); saveat=0.1, save_idxs=1)
        period, amplitude = getPerAmp(sol)

        push!(extended_dataset, (ind=ind, fit=fitness, per=period, amp=amplitude))
    end
    extended_dataset_df = DataFrame(extended_dataset)
    extended_dataset_df = df[df.per .!= 0.0, :] # remove individuals with period of 0
    extended_dataset_df.freq = 100 ./ df.per # convert period column to frequency
    return extended_dataset_df
end

df = analyze_population(filtered_fitpops)



#! SCATTER PLOTS
begin
    n_params = length(param_names)

    # Create the plots
    ka_params = ["ka1", "ka2", "ka3", "ka4", "ka7"]
    kb_params = ["kb1", "kb2", "kb3", "kb4", "kb7"]
    kcat_params = ["kcat1", "kcat7"]

    # Helper function to get the index of the parameter in the ind column
    get_param_index = (param) -> findfirst(==(param), param_names)

    # Plot ka parameters
    ka_plot = plot(layout = (2, 5), size = (1000, 400), dpi=200, legend=false)
    for (i, param) in enumerate(ka_params)
        scatter!(ka_plot, [x[get_param_index(param)] for x in df.ind], df[:, :amp], subplot = i, xlabel = param, ylabel = "Amplitude (uM)", color = :red, markersize=2)
        scatter!(ka_plot, [x[get_param_index(param)] for x in df.ind], df[:, :per], subplot = i + 5, xlabel = param, ylabel = "Period (s)", color = :blue, markersize=2)
    end

    # Plot kb parameters
    kb_plot = plot(layout = (2, 5), size = (1000, 400), dpi=200, legend=false)
    for (i, param) in enumerate(kb_params)
        scatter!(kb_plot, [x[get_param_index(param)] for x in df.ind], df[:, :amp], subplot = i, xlabel = param, ylabel = "Amplitude (uM)", color = :red, markersize=2)
        scatter!(kb_plot, [x[get_param_index(param)] for x in df.ind], df[:, :per], subplot = i + 5, xlabel = param, ylabel = "Period (s)", color = :blue, markersize=2)
    end

    # Plot kcat parameters
    kcat_plot = plot(layout = (2, 2), size = (800, 400), dpi=200, legend=false)
    for (i, param) in enumerate(kcat_params)
        scatter!(kcat_plot, [x[get_param_index(param)] for x in df.ind], df[:, :amp], subplot = i, xlabel = param, ylabel = "Amplitude (uM)", color = :red, markersize=2)
        scatter!(kcat_plot, [x[get_param_index(param)] for x in df.ind], df[:, :per], subplot = i + 2, xlabel = param, ylabel = "Period (s)", color = :blue, markersize=2)
    end

    # Plot y parameter
    y_plot = plot(layout = (2, 1), size = (800, 400), dpi=200, legend=false)
    y_index = get_param_index("y")
    scatter!(y_plot, [x[y_index] for x in df.ind], df[:, :amp], subplot = 1, xlabel = "y (V/A)", ylabel = "Amplitude (uM)", color = :red, markersize=2)
    scatter!(y_plot, [x[y_index] for x in df.ind], df[:, :per], subplot = 2, xlabel = "y (V/A)", ylabel = "Period (s)", color = :blue, markersize=2)

    # Display the plots or save them using savefig("filename.png")
    display(ka_plot); savefig(ka_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/scatter_ka.png")
    display(kb_plot); #savefig(kb_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/scatter_kb.png")
    display(kcat_plot); #savefig(kcat_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/scatter_kcat.png")
    display(y_plot); #savefig(y_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/scatter_y.png")
end


#! Pearson correlation heatmap
begin
    # Extract the parameter values from the DataFrame
    params = hcat([[x[i] for x in df.ind] for i in 1:n_params]...)

    # Calculate the correlation matrix
    cor_matrix = cor(params, dims=1)

    # Create a heatmap of the correlation matrix
    heatmap_plot = heatmap(param_names, param_names, cor_matrix, aspect_ratio=1, clims=(-1, 1), color=:coolwarm)
    title!("Parameter Correlations")
    xlabel!("Parameters")
    ylabel!("Parameters")

    # Display the heatmap or save it using savefig("filename.png")
    display(heatmap_plot)
    savefig(heatmap_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/correlation_heatmap.png")
end


#! SURFACE PLOTS
function plot_parameter_surface(df::DataFrame, param1::String, param2::String, n_points::Int=100)
    # Get the indices of the input parameters
    param1_index = findfirst(==(param1), param_names)
    param2_index = findfirst(==(param2), param_names)

    # Calculate the parameter value arrays based on the DataFrame
    param1_values = [x[param1_index] for x in df.ind]
    param2_values = [x[param2_index] for x in df.ind]

    # Calculate the mean and standard deviation for each parameter
    param1_mean, param1_std = mean(param1_values), std(param1_values)
    param2_mean, param2_std = mean(param2_values), std(param2_values)

    # Calculate the parameter ranges
    param1_range = range(param1_mean - param1_std, param1_mean + param1_std, length=n_points)
    param2_range = range(param2_mean - param2_std, param2_mean + param2_std, length=n_points)

    # Preallocate matrices for amplitude and period
    amplitude = zeros(n_points, n_points)
    period = zeros(n_points, n_points)

    # Loop over the parameter ranges and evaluate the model
    for (i, p1) in enumerate(param1_range), (j, p2) in enumerate(param2_range)
        # Update the model parameters
        new_params = deepcopy(df.ind[1])
        new_params[param1_index] = p1
        new_params[param2_index] = p2

        # Solve the model and calculate amplitude and period
        sol = solve(remake(prob, p=new_params), save_idxs=1, verbose=false) # Replace with your model-solving function
        per, amp = getPerAmp(sol)

        # Store the results in the matrices
        amplitude[i, j] = amp
        period[i, j] = per
    end

    # Create the surface plots
    p1 = surface(param1_range, param2_range, amplitude, xlabel=param1, ylabel=param2, zlabel="Amplitude (uM)", alpha=0.7)
    p2 = surface(param1_range, param2_range, period, xlabel=param1, ylabel=param2, zlabel="Period (s)", alpha=0.7)

    # Combine the plots into side-by-side subplots
    combined_plot = plot(p1, p2, layout=(1, 2), size=(1000, 500), dpi=200)

    # Display the plots or save them using savefig("filename.png")
    display(combined_plot)
    savefig(combined_plot, "OscillatorPaper/MaggieModelReduction/FIGURES/SurfacePlots/$(param1)_$(param2).png")
end

# Call the function with the desired parameters and DataFrame
plot_parameter_surface(df, "ka3", "y")




#! Dimensional reduction 
# Preprocess the data
X = hcat(df.ind...)'
X_centered = X .- mean(X, dims=1)

# PCA
pca = fit(PCA, X_centered', maxoutdim=2)
X_pca = MultivariateStats.transform(pca, X_centered')
X_pca = X_pca'
# t-SNE
X_tsne = tsne(X, 2)

# UMAP
umap_model = UMAP_(X', 2, n_neighbors=14)
X_umap = umap_model.embedding'

# Visualization
p1 = scatter(X_pca[:, 1], X_pca[:, 2], zcolor=df.amp, title="PCA", legend=false, colorbar=true)
p2 = scatter(X_tsne[:, 1], X_tsne[:, 2], zcolor=df.amp, title="t-SNE", legend=false, colorbar=true)
p3 = scatter(X_umap[:, 1], X_umap[:, 2], zcolor=df.amp, title="UMAP", legend=false, colorbar=true)
plot(p1, p2, p3, layout=(3, 1), size=(900, 1100))


function plot_reduced_spaces(df::DataFrame, X_pca::AbstractArray, X_tsne::AbstractArray, X_umap::AbstractArray)
    # Normalize amplitude and period for color mapping
    amplitude = df[:, :amp] ./ maximum(df[:, :amp])
    period = df[:, :per] ./ maximum(df[:, :per])

    plots_amp = []
    plots_per = []

    methods = ["PCA", "t-SNE", "UMAP"]
    reduced_spaces = [X_pca, X_tsne, X_umap]

    for (i, method) in enumerate(methods)
        X = reduced_spaces[i]
        if size(X, 2) < 2
            println("Warning: $method has fewer than 2 dimensions. Skipping this method.")
            continue
        end

        scatter_amp = scatter(X[:, 1], X[:, 2], zcolor=amplitude, xlabel="$(method) 1", ylabel="$(method) 2", title=method, legend=:none, colorbar=true)
        scatter_per = scatter(X[:, 1], X[:, 2], zcolor=period, xlabel="$(method) 1", ylabel="$(method) 2", title=method, legend=:none, colorbar=true)

        push!(plots_amp, scatter_amp)
        push!(plots_per, scatter_per)
    end

    # Combine the plots into side-by-side subplots with amplitude and period as color
    combined_plot_amp = plot(plots_amp..., title="Amplitude (uM)", layout=(length(plots_amp),1), size=(600,1000), dpi=200)
    combined_plot_per = plot(plots_per..., title= "Period (s)", layout=(length(plots_per),1), size=(600,1000), dpi=200)

    display(combined_plot_amp); savefig(combined_plot_amp, "OscillatorPaper/MaggieModelReduction/FIGURES/DimensionReduction/Amplitude.png")
    display(combined_plot_per); savefig(combined_plot_per, "OscillatorPaper/MaggieModelReduction/FIGURES/DimensionReduction/Period.png")
end

plot_reduced_spaces(df, X_pca, X_tsne, X_umap)




#! Sensitivity Analysis
function create_sensitivity_analysis_func(prob, eval_func)
    function inner_sensitivity_analysis(P)
        prob_func(prob, i, repeat) = remake(prob; p=P[:, i])
        ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
        ntraj = size(P, 2)
        ensol = solve(ensemble_prob, EnsembleThreads(); saveat=0.1, save_idxs = 1, maxiters=10000, verbose=false, trajectories=ntraj)

        out = zeros(ntraj)
        for i in 1:ntraj
            out[i] = eval_func(ensol[i])
        end
        out
    end
    return inner_sensitivity_analysis
end

peramp_sensitivity_analysis_func = create_sensitivity_analysis_func(prob, getPerAmp)
oscillation_sensitivity_analysis_func = create_sensitivity_analysis_func(prob, CostFunction)

function extract_param_bounds(param_values::Dict{String, Dict{String, Float64}})
    lb = Float64[]
    ub = Float64[]
    for (key, value) in param_values
        push!(lb, value["min"])
        push!(ub, value["max"])
    end
    return lb, ub
end


samples = 1000
sampler = SobolSample()
lb, ub = extract_param_bounds(param_values)
A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)

sobol_result = gsa(oscillation_sensitivity_analysis_func, Sobol(), A, B, batch=true)
first_order_indices = sobol_result.S1
total_order_indices = sobol_result.ST

function plot_sobol_indices(first_order_indices, total_order_indices, param_names)
    p1 = bar([first_order_indices], xlabel="Parameters", ylabel="First-order Sobol indices", legend=:none)
    xticks!(p1, (collect(1:length(param_names)), param_names))

    p2 = bar([total_order_indices], xlabel="Parameters", ylabel="Total-order Sobol indices", legend=:none)
    xticks!(p2, (collect(1:length(param_names)), param_names))

    sobolplot = plot(p1, p2, layout=(2, 1), size=(1000, 600))
    display(sobolplot); savefig(sobolplot, "OscillatorPaper/MaggieModelReduction/FIGURES/SensitivityAnalysis/SobolIndices.png")
end

plot_sobol_indices(first_order_indices, total_order_indices, param_names)




## ! BIFURCATION ANALYSIS ##
# * define bifurcation parameter and interval to vary it over
# TODO: change plot var to amplitude or period
bif_par = :ka3           # bifurcation parameter
p_span = (0.1, 20.)    # interval to vary ka3 over
plot_var = :L          # we will plot L vs ka3
p_bstart = copy(psym)
p_bstart[findfirst(x -> x.first == bif_par, psym)] = :ka3 => p_span[1]

# * define the ODEProblem with derivative and Jacobian functions
oprob = ODEProblem(fullrn, usym, (0.0,0.0), p_bstart; jac=true) # ? don't know why tspan is (0.0,0.0)? 
F = (u,p) -> oprob.f(u, p, 0) # ? don't know why t is 0?
J = (u,p) -> oprob.f.jac(u, p, 0)

# * get ka3 and L as symbolic variables
@unpack ka3, L = fullrn

# * find their indices in oprob.p and oprob.u0 respectively
bif_idx  = findfirst(isequal(ka3), parameters(fullrn))
plot_idx = findfirst(isequal(L), species(fullrn))

u00 = [0.0, 3.0, 0.5, 0.3, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# * define the Bifurcation problem
bprob1 = BifurcationProblem(F, u00, oprob.p, (@lens _[bif_idx]);
                           recordFromSolution = (x, p) -> x[plot_idx], J = J)

# * define the options for the continuation method
bopts = ContinuationPar(dsmax = 0.01,          # Max arclength in PACM.
        dsmin = 1e-4,          # Min arclength in PACM.
        ds = 0.002,            # Initial (positive) arclength in PACM.
        maxSteps = 100000,     # Max number of steps.
        pMin = p_span[1],      # Min p-val (if hit, the method stops).
        pMax = p_span[2],      # Max p-val (if hit, the method stops).
        detectBifurcation = 3,
        nInversion = 6, 
        maxBisectionSteps = 25,
        nev = 4,) # Value in {0,1,2,3}

# * We compute periodic orbits of (E) using the Trapezoid method. We start with two helper functions that record and plot the periodic orbits.
argspo = (recordFromSolution = (x, p) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p)))
	end,
	plotSolution = (x, p; k...) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
		# plot!(br; subplot = 1, putspecialptlegend = false)
	end)

# * this is the function which builds probtrap from sol
probtrap, ci = BK.generateCIProblem(PeriodicOrbitTrapProblem(M = 151;
	jacobian = :DenseAD, updateSectionEveryStep = 0),
	bprob1, sol, 2.)

opts_po_cont = setproperties(opts_br, maxSteps = 50, tolStability = 1e-8)
brpo_fold = continuation(probtrap, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)

scene = plot(brpo_fold)



# * trying parallel standard shooting
record_sh = recordFromSolution = (x, p) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, set(par_pop, p.prob.lens, p.p))
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, set(par_pop, p.prob.lens, p.p)))
	end
plot_sh  = (x, p; k...) -> begin
	xtt = BK.getPeriodicOrbit(p.prob, x, set(par_pop, p.prob.lens, p.p))
	plot!(xtt.t, xtt[1,:]; label = "x", k...)
	plot!(xtt.t, xtt[2,:]; label = "y", k...)
	# plot!(br; subplot = 1, putspecialptlegend = false)
	end

probsh, cish = generateCIProblem( ShootingProblem(M=3),
	bprob, prob, sol, 2.; alg = Rodas5())

opts_po_cont = setproperties(opts_br, maxSteps = 50, saveEigenvectors = true, detectLoop = true, tolStability = 1e-3)
br_fold_sh = continuation(probsh, cish, PALC(tangent = Bordered()), opts_po_cont;
	verbosity = 3, plot = true,
	recordFromSolution = record_sh,
	plotSolution = plot_sh,)

scene = plot(br_fold_sh)


# * perform bifurcation analysis
# TODO: need to find right continuation algorithm
bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> bopts)

plot(bf; xlabel = string(bif_par), ylabel = string(plot_var))






#! Plot oscillatory regions of 2D parameter space with contour or heatplot 
param_ranges = ("kb7" => range(start = 0.01, stop = 100.0; length = 300), "V/A" => range(start = 100., stop = 10000.0, length = 300))

#* Helper function to assign units based on string
function assign_units(param_name::String)
    if occursin("ka",param_name)
        return u"1/(µM*s)"
    elseif occursin("kb",param_name)
        return u"1/s"
    elseif occursin("kcat",param_name)
        return u"1/s"
    elseif param_name == "V/A"
        return u"nm"
    end
end

#* Define the function to evaluate the 2D parameter space for oscillations
function evaluate_2D_parameter_space(param_ranges::Tuple, prob::ODEProblem)
    p1name, p2name = param_ranges[1].first, param_ranges[2].first

    #? Get indices of parameters in the parameter vector
    p1idx = findfirst(x -> x == p1name, param_names)
    p2idx = findfirst(x -> x == p2name, param_names)

    #? Get units of parameters
    p1_units = assign_units(p1name)
    p2_units = assign_units(p2name)
    
    #? Get evaluation function
    evalfunc = make_fitness_function(eval_fitness_catcherrors, prob)
    
    #? Create a copy of the parameter vector
    newparams = copy(prob.p)
    
    #? Define a function to pass modified parameters to the evaluation function
    function is_oscillatory(newvals)
        newparams[p1idx] = newvals[1]
        newparams[p2idx] = newvals[2]
    
        evalfunc(newparams) #< -0.5 ? true : false
    end

    #? Get ranges of parameters
    p1_range = param_ranges[1].second
    p2_range = param_ranges[2].second

    #? Get dimensions of parameter space
    n_p1 = length(p1_range)
    n_p2 = length(p2_range)

    #? Create arrays to store results
    oscillatory_points = Array{Float64}(undef, n_p1, n_p2)
    p1_points = Vector{Float64}(undef, n_p1)
    p2_points = Vector{Float64}(undef, n_p2)

    #? Progress bar
    progress = Progress(n_p1*n_p2, dt = 0.1)

    #? Evaluate parameter space
    for (i, p1_val) in enumerate(p1_range)
        p1_points[i] = p1_val
        for (j, p2_val) in enumerate(p2_range)
            p2_points[j] = p2_val
            oscillatory_points[i, j] = is_oscillatory((p1_val, p2_val))
            next!(progress)
        end
    end

    return p1_points*p1_units, p2_points*p2_units, -oscillatory_points
end


function plot_oscillatory_regions(p1_points, p2_points, oscillatory_points, param_ranges)
    p1name, p2name = param_ranges[1].first, param_ranges[2].first
    p = contour(p1_points, p2_points, oscillatory_points', xlabel=p1name, ylabel=p2name, 
                title="Oscillatory Regions of Parameter Space", color=:vik, legend=:none, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    display(p)
    if any(i -> i == "V/A", (p1name,p2name))
        p1name == "V/A" ? p1name = "VA" : p2name = "VA"
    end
    savefig(p, "OscillatorPaper/MaggieModelReduction/FIGURES/ParameterSpace/OscillatoryRegions_$(p1name)_$(p2name).png")
end

# Usage example
p1_points, p2_points, oscillatory_points = evaluate_2D_parameter_space(param_ranges, prob)
plot_oscillatory_regions(p1_points, p2_points, oscillatory_points, param_ranges)

testp = copy(prob.p)
testp[6] = 5.0
testp[end] = 3000.0
testsol = solve(remake(prob, p = testp))
plot(testsol)




#! Reachability analysis with two fixed parameters. Measure volume of oscillatory region.
# The optimization function using a genetic algorithm (returns a filtered list of oscillatory points)
function optimize_parameters(param_values::OrderedDict, new_fitness_function::Function)

    # Run your genetic algorithm optimization code with the updated param_values
    begin
        population_size = 8000
        pop = generate_population(param_values, population_size)
        constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
        opts = Evolutionary.Options(abstol=1e-12, reltol=1e-10, successive_f_tol = 2, iterations=6, store_trace = true, 
                show_trace=true, show_every=1, parallelization=:thread)
        common_range = 0.5; valrange = fill(common_range, length(param_values))
        mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
                crossover = TPX, crossoverRate = 0.5,
                mutation  = BGA(valrange, 2), mutationRate = 0.7)

        #? Convergence progress bar
        # progress = ProgressThresh(opts.abstol, "Converging: ")

        #* Optimization
        result = Evolutionary.optimize(new_fitness_function, constraints, mthd, pop, opts)
    end
    fitpops = [gen.metadata["populationmap"] for gen in result.trace]
    fitpops = reduce(vcat, fitpops)
    
    # Apply fitness threshold
    filter!(x -> x.fit <= -0.4, fitpops)
    # Then extract the oscillatory points as you mentioned:
    oscillatory_points = [x.ind for x in fitpops]
    
    return oscillatory_points
end

"""Make new fitness function that includes the fixed parameters"""
function make_new_fitness_function(func::Function, prob::ODEProblem, fixed_params::OrderedDict)
    function new_fitness_function(tracker::PeriodAmplitudes, p::Vector{Float64})
        p_combined = Vector{Float64}(undef, length(p) + length(fixed_params))
        i_fixed, i_var = 1, 1

        for (i, key) in enumerate(keys(param_values))
            if key in keys(fixed_params)
                p_combined[i] = fixed_params[key]
                i_fixed += 1
            else
                p_combined[i] = p[i_var]
                i_var += 1
            end
        end

        return func(tracker, p_combined, prob)
    end
    return new_fitness_function
end

"""Monte carlo sampling to estimate bounding volume"""
function monte_carlo_volume(points::Vector{Vector{T}}, param_values::OrderedDict, n_samples::Int) where {T<:Real}
    n_points, dim = length(points), length(points[1])

    # Calculate the total volume using the parameter constraint ranges
    total_volume = prod([param_values[key]["max"] - param_values[key]["min"] for key in keys(param_values)])
    @info "Total parameter volume: $total_volume"

    # Generate random points within the bounds specified by param_values
    random_points = [rand() * (param_values[key]["max"] - param_values[key]["min"]) + param_values[key]["min"] for key in keys(param_values) for _ in 1:n_samples]

    # Reshape the random_points into an array of vectors
    random_points = reshape(random_points, (dim, n_samples))

    # Calculate distances between random points and input points
    distances = [norm(points[j] .- random_points[:, i]) for i in 1:n_samples, j in 1:n_points]

    # Find the minimum distance from each random point to any of the input points
    min_distances = minimum(distances, dims=2)

    # Check if each random point is within the maximum distance between the input points
    in_hull = sum(min_distances .<= maximum([norm(points[i] .- points[j]) for i in 1:n_points for j in 1:n_points]) * 2)

    estimated_volume = total_volume * (in_hull / n_samples)

    return estimated_volume
end


monte_carlo_volume(fitinds, param_values, 10000)

function generate_random_points(param_values, n_samples)
    random_points = Vector{Vector{Float64}}(undef, 0)
    for _ in 1:n_samples
        point = []
        for key in keys(param_values)
            min_val = param_values[key]["min"]
            max_val = param_values[key]["max"]
            element = rand(min_val:0.001:max_val)
            push!(point, element)
        end
        push!(random_points, point)
    end
    return random_points
end

random_points = generate_random_points(param_values, 1000)


function reachability_analysis(param_values::OrderedDict{String, Dict{String, Float64}}, fixed_param_pairs::Base.Generator, prob::ODEProblem)
    # results = OrderedDict{Vector{String}, Tuple{Float64, Int, Float64, Float64, Float64}}()
    volumes = []
    num_oscillatory_points = []
    avg_fitnesses = []
    avg_periods = []
    avg_amplitudes = []

    loopprogress = Progress(length(collect(fixed_param_pairs)), desc ="Looping thru fixed pairs: " , color=:green)

    @info "Starting loop"

    for param_pair in fixed_param_pairs
        @info "Fixed parameters pair: $param_pair"

        fixed_params = OrderedDict(param_pair .=> [fixed_param_values[p] for p in param_pair])
        variable_param_values = deepcopy(param_values)
        delete!.(Ref(variable_param_values), keys(fixed_params))

        # Create a new instance of PeriodAmplitudes for each param_pair loop
        tracker = PeriodAmplitudes()

        new_fitness_function = make_new_fitness_function(eval_fitness_catcherrors, prob, fixed_params)

        # Wrap the fitness_function call to work with the Evolutionary.jl optimize function
        wrapped_fitness_function(params) = new_fitness_function(tracker, params)

        # Run the optimization function to get the oscillatory points
        oscillatory_points = optimize_parameters(variable_param_values, wrapped_fitness_function)::Vector{Vector{Float64}} #? Vector of oscillatory points

        # Compute convex hull volume of the oscillatory region in parameter space
        volume = LazySets.volume(oscillatory_points)
        variable_param_keys = collect(keys(variable_param_values))
        min_values = OrderedDict(p => Inf for p in variable_param_keys)
        max_values = OrderedDict(p => -Inf for p in variable_param_keys)

        for point in oscillatory_points
            vals = point.ind
            for (p, val) in zip(variable_param_keys, vals)
                min_values[p] = min(min_values[p], val)
                max_values[p] = max(max_values[p], val)
            end
        end

        # Compute the volume of the oscillatory region in parameter space
        @show normalized_volume = prod((max(0, max_values[p] - min_values[p]) / (param_values[p]["max"] - param_values[p]["min"])) for p in variable_param_keys)
        if normalized_volume == 0.0
            @warn "Oscillatory region has zero volume"
            continue
        end

        # Calculate the number of oscillatory points
        num_points = length(oscillatory_points)
        @info "Number of oscillatory points: $num_points"

        # Calculate the average fitness value
        avg_fitness = mean(point.fit for point in oscillatory_points)

        # Calculate the average period and amplitude
        # avg_period = mean(haskey(tracker.peramps, point.ind) ? tracker.peramps[point.ind][1] : 0 for point in oscillatory_points)
        avg_period = mean(x -> x[2][1], tracker.peramps)
        @info "Average period: $avg_period"
        # avg_amplitude = mean(haskey(tracker.peramps, point.ind) ? tracker.peramps[point.ind][2] : 0 for point in oscillatory_points)
        avg_amplitude = mean(x -> x[2][2], tracker.peramps)
        @info "Average amplitude: $avg_amplitude"


        # Store the results for the given fixed parameter combination
        # results[param_pair] = (normalized_volume, num_points, avg_fitness, avg_period, avg_amplitude)
        push!(volumes, normalized_volume)
        push!(num_oscillatory_points, num_points)
        push!(avg_fitnesses, avg_fitness)
        push!(avg_periods, avg_period)
        push!(avg_amplitudes, avg_amplitude)


        next!(loopprogress)
    end
    # Convert results to DataFrame
    results_df = DataFrame(combo = collect(fixed_param_pairs), volume = volumes, num_points = num_oscillatory_points, avg_fitness = avg_fitnesses, avg_period = avg_periods, avg_amplitude = avg_amplitudes)
    return results_df
end

function visualize_reachability(results::Dict{Vector{String}, Tuple{Float64, Int, Float64, Float64, Float64}})
    x_vals = Float64[]
    y_vals = Float64[]
    x_labels = String[]
    y_labels = String[]
    reachability_vals = Float64[]
    num_points_vals = Float64[]
    avg_fitness_vals = Float64[]
    avg_period_vals = Float64[]
    avg_amplitude_vals = Float64[]

    for (key, val) in results
        push!(x_vals, fixed_param_values[key[1]])
        push!(y_vals, fixed_param_values[key[2]])
        push!(x_labels, key[1])
        push!(y_labels, key[2])
        push!(reachability_vals, val[1])
        push!(num_points_vals, val[2])
        push!(avg_fitness_vals, val[3])
        push!(avg_period_vals, val[4])
        push!(avg_amplitude_vals, val[5])
    end

    p1 = scatter(x_vals, y_vals, reachability_vals, xlabel=x_labels[1], ylabel=y_labels[1], zlabel="Volume", title="Reachability Volume", color=:viridis, legend=:none)
    p2 = scatter(x_vals, y_vals, num_points_vals, xlabel=x_labels[1], ylabel=y_labels[1], zlabel="Number of Points", title="Number of Oscillatory Points", color=:viridis, legend=:none)
    p3 = scatter(x_vals, y_vals, avg_fitness_vals, xlabel=x_labels[1], ylabel=y_labels[1], zlabel="Average Fitness", title="Average Fitness Value", color=:viridis, legend=:none)
    p4 = scatter(x_vals, y_vals, avg_period_vals, xlabel=x_labels[1], ylabel=y_labels[1], zlabel="Average Period", title="Average Period", color=:viridis, legend=:none)
    p5 = scatter(x_vals, y_vals, avg_amplitude_vals, xlabel=x_labels[1], ylabel=y_labels[1], zlabel="Average Amplitude", title="Average Amplitude", color=:viridis, legend=:none)
    
    plot(p1, p2, p3, p4, p5, layout=(1, 5), size=(2000, 400))
end




# function visualize_reachability(results::Dict{Vector{String}, Tuple{Float64, Int, Float64}})
#     x_vals = Float64[]
#     y_vals = Float64[]
#     reachability_vals = Float64[]
#     num_points_vals = Float64[]
#     avg_fitness_vals = Float64[]

#     for (key, val) in results
#         push!(x_vals, fixed_param_values[key[1]])
#         push!(y_vals, fixed_param_values[key[2]])
#         push!(reachability_vals, val[1])
#         push!(num_points_vals, val[2])
#         push!(avg_fitness_vals, val[3])
#     end

#     p1 = scatter(x_vals, y_vals, reachability_vals, xlabel=key[1], ylabel=key[2], zlabel="Volume", title="Reachability Volume", color=:viridis, legend=:none)
#     p2 = scatter(x_vals, y_vals, num_points_vals, xlabel=key[1], ylabel=key[2], zlabel="Number of Points", title="Number of Oscillatory Points", color=:viridis, legend=:none)
#     p3 = scatter(x_vals, y_vals, avg_fitness_vals, xlabel=key[1], ylabel=key[2], zlabel="Average Fitness", title="Average Fitness Value", color=:viridis, legend=:none)
    
#     plot(p1, p2, p3, layout=(1, 3), size=(1200, 400))
# end



# Perform the reachability analysis and visualize the results
fixed_param_values = OrderedDict(param_names[i] => p[i] for i in eachindex(param_names))
fixed_param_names = collect(keys(fixed_param_values))
fixed_param_pairs = combinations(fixed_param_names, 2)

#! Run the reachability analysis
results = reachability_analysis(param_values, fixed_param_pairs, prob)


#! Visualize the results
visualize_reachability(results)


function heatmap_normalized_volumes(results)
    x_labels = sort(unique([key[1] for key in keys(results)]))
    @info "x_labels: $x_labels"
    y_labels = sort(unique([key[2] for key in keys(results)]))
    @info "y_labels: $y_labels"

    num_x_labels = length(x_labels)
    num_y_labels = length(y_labels)

    x_indices = Dict(x_labels[i] => i for i in 1:num_x_labels)
    @info "x_indices: $x_indices"
    y_indices = Dict(y_labels[i] => i for i in 1:num_y_labels)
    @info "y_indices: $y_indices"

    volumes = zeros(num_x_labels, num_y_labels)

    for (key, val) in results
        x_index = x_indices[key[1]]
        y_index = y_indices[key[2]]
        volumes[x_index, y_index] = val[1]
    end

    @show min_volume = maximum(volumes)
    normalized_volumes = volumes ./ min_volume
    clamped_normalized_volumes = clamp.(normalized_volumes .- 1, -1 + eps(Float64), Inf)
    log_normalized_volumes = log1p.(clamped_normalized_volumes)

    heatmap(x_labels, y_labels, log_normalized_volumes, xlabel="Fixed Parameter 1", ylabel="Fixed Parameter 2", title="Log Normalized Reachability Volume")
end



heatmap_normalized_volumes(results)



function scatter_volume_vs_points(results)
    volumes = [val[1] for val in values(results)]
    num_points = [val[2] for val in values(results)]
    
    scatter(volumes, num_points, xlabel="Reachability Volume", ylabel="Number of Oscillatory Points", title="Volume vs Oscillatory Points")
end

scatter_volume_vs_points(results)

function barplot_average_fitness(results)
    param_pairs = [string(key[1], ", ", key[2]) for key in keys(results)]
    avg_fitness_values = [val[3] for val in values(results)]
    
    bar(param_pairs, avg_fitness_values, xlabel="Fixed Parameter Pairs", ylabel="Average Fitness", title="Average Fitness for Fixed Parameter Pairs")
end

barplot_average_fitness(results)

function bubble_plot_reachability_volume(results)
    x_vals = [fixed_param_values[key[1]] for key in keys(results)]
    y_vals = [fixed_param_values[key[2]] for key in keys(results)]
    reachability_volumes = [val[1] for val in values(results)]
    
    # Normalize marker sizes to avoid extremely small bubbles
    normalized_sizes = (reachability_volumes .- minimum(reachability_volumes)) ./ (maximum(reachability_volumes) - minimum(reachability_volumes))
    markersizes = 5 .+ normalized_sizes .* 50

    scatter(x_vals, y_vals, markersize=markersizes, xlabel="Fixed Parameter 1", ylabel="Fixed Parameter 2", title="Reachability Volume (Bubble Plot)", legend=:none)
end

bubble_plot_reachability_volume(results)

function matrix_plot_period_amplitude(results)
    unique_x_labels = unique([key[1] for key in keys(results)])
    unique_y_labels = unique([key[2] for key in keys(results)])

    matrix_period = zeros(length(unique_x_labels), length(unique_y_labels))
    matrix_amplitude = zeros(length(unique_x_labels), length(unique_y_labels))

    for (key, val) in results
        i = findfirst(isequal(key[1]), unique_x_labels)
        j = findfirst(isequal(key[2]), unique_y_labels)
        matrix_period[i, j] = val[4]
        matrix_amplitude[i, j] = val[5]
    end

    p1 = heatmap(unique_x_labels, unique_y_labels, matrix_period, xlabel="Fixed Parameter 1", ylabel="Fixed Parameter 2", title="Average Period")
    p2 = heatmap(unique_x_labels, unique_y_labels, matrix_amplitude, xlabel="Fixed Parameter 1", ylabel="Fixed Parameter 2", title="Average Amplitude")

    plot(p1, p2, layout=(1, 2), size=(1000, 400))
end

matrix_plot_period_amplitude(results)

function scatter_fitness_vs_volume(results)
    avg_fitness_values = [val[3] for val in values(results)]
    volumes = [val[1] for val in values(results)]

    scatter(avg_fitness_values, volumes, xlabel="Average Fitness", ylabel="Reachability Volume", title="Fitness vs Reachability Volume")
end

scatter_fitness_vs_volume(results)

function scatter_period_vs_amplitude(results)
    avg_period_values = [val[4] for val in values(results)]
    avg_amplitude_values = [val[5] for val in values(results)]

    scatter(avg_period_values, avg_amplitude_values, xlabel="Average Period", ylabel="Average Amplitude", title="Period vs Amplitude")
end

scatter_period_vs_amplitude(results)

function stacked_barplot_period_amplitude(results)
    param_pairs = [string(key[1], ", ", key[2]) for key in keys(results)]
    avg_period_values = [val[4] for val in values(results)]
    avg_amplitude_values = [val[5] for val in values(results)]

    bar(param_pairs, [avg_period_values, avg_amplitude_values], label=["Average Period" "Average Amplitude"], xlabel="Fixed Parameter Pairs", ylabel="Value", title="Average Period and Amplitude for Fixed Parameter Pairs")
end

stacked_barplot_period_amplitude(results)

# using Plots
pyplot()

function radar_plot_all_metrics(results)
    param_pairs = [string(key[1], ", ", key[2]) for key in keys(results)]
    normalized_volumes = [val[1] for val in values(results)]
    num_points = [val[2] for val in values(results)]
    avg_fitness = [val[3] for val in values(results)]
    avg_period = [val[4] for val in values(results)]
    avg_amplitude = [val[5] for val in values(results)]

    data = hcat(normalized_volumes, num_points, avg_fitness, avg_period, avg_amplitude)
    labels = ["Normalized Volume", "Num Points", "Avg Fitness", "Avg Period", "Avg Amplitude"]

    @df DataFrame(data, Symbol.(labels)) radarplot(labels, data, title="Radar Plot for All Metrics", legend=:topleft)
end

radar_plot_all_metrics(results)