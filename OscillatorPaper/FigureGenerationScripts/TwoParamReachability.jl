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
    # using Unitful
    # using Unitful: µM, nm, s
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
    # plotlyjs()
    # gr()
    # using Base.Threads


    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")
end


#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
    #? Parameter list
    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #? Initial condition list
    usym = [:L => 0.0, :K => 0.5, :P => 0.3, :A => 2.0, :Lp => 3.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    fullprob = ODEProblem(fullrn, u0, tspan, p)
    sol = solve(fullprob, saveat=0.1, save_idxs=1)

    #? Plot the results
    # plot(sol)
end


"""Returns the function factory for the cost function, referencing the ODE problem, tracker, and fixed inputs with closure"""
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_pair::Vector{ConstraintRange}, pair_idxs::Tuple{Int,Int})
    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_pair))

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if fixed_inputs_inserted < length(fixed_input_pair) && i == pair_idxs[fixed_inputs_inserted + 1]
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_pair[fixed_inputs_inserted + 1].nominal
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[i - fixed_inputs_inserted]
            end
        end

        return evalfunc(new_input, prob)
    end
    return fitness_function
end


"""Monte carlo sampling to estimate bounding volume"""
function monte_carlo_volume(points::Vector{Vector{Float64}}, constraint::ConstraintType, n_samples::Int = 10000)
    n_points, dim = length(points), length(points[1])

    # Calculate the total volume using the parameter constraint ranges
    total_volume = prod([conrange.max - conrange.min for conrange in constraint.ranges])
    @info "Total solution volume: $total_volume"

    # Generate random points within the bounds specified by param_values
    random_points = [rand() * (conrange.max - conrange.min) + conrange.min for conrange in constraint.ranges for _ in 1:n_samples]

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


"""Compare the allowed solution space for when each pair of inputs is fixed"""
function reachability_analysis(constraints::ConstraintType, prob::ODEProblem)

    # Get the nominal values of the inputs depending on whether they are parameters or initial conditions. 
    # nominal_values = constraints isa ParameterConstraints ? SVector{13}(prob.p) : nominal_values = SVector{4}(prob.u0[1:4]) 

    #Now, map each input name to its nominal value within fixed_input_combos. This is what we will loop through
    fixed_pair_combos = combinations(constraints.ranges, 2) # Get all combinations of name:nominalvalue.

    loopprogress = Progress(length(collect(fixed_pair_combos)), desc ="Looping thru fixed pairs: " , color=:green)


    #* Make a results dictionary where fixedpair => (volume, avg_period, avg_amplitude, num_oscillatory_points)
    results = Dict{Tuple{Symbol, Symbol}, NamedTuple{(:volume, :avg_period, :avg_amplitude, :num_oscillatory_points), Tuple{Float64, Float64, Float64, Int}}}() 

    for fixedpair in fixed_pair_combos
        @info "Fixed input pair: $(fixedpair[1].name), $(fixedpair[2].name)"

        # Make a copy of the constraints to modify
        variable_constraints = deepcopy(constraints)


        # Remove the fixed parameters from the variable constraints
        filter!(x -> x.name != fixedpair[1].name && x.name != fixedpair[2].name, variable_constraints.ranges)

        fixed_ga_problem = GAProblem(variable_constraints, prob)

        fixedpair_idxs = find_indices(fixedpair, constraints.ranges) # Get the indices of the fixed inputs.

        #* Create a closure for the fitness function that includes the fixed inputs
        make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixedpair, fixedpair_idxs)

        # Run the optimization function to get the oscillatory points
        oscillatory_points = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 8) #? Vector of oscillatory points
        # return oscillatory_points
        #! START HERE ##
        # Compute convex hull volume of the oscillatory region in parameter space
        volume = monte_carlo_volume([x.ind for x in oscillatory_points], constraints)
 

        # Compute the volume of the oscillatory region in parameter space
        # @show normalized_volume = prod((max(0, max_values[p] - min_values[p]) / (param_values[p]["max"] - param_values[p]["min"])) for p in variable_param_keys)
        # if normalized_volume == 0.0
        #     @warn "Oscillatory region has zero volume"
        #     continue
        # end

        # Calculate the number of oscillatory points
        num_points = length(oscillatory_points)
        @info "Number of oscillatory points: $num_points"

        # Calculate the average period and amplitude
        avg_period = mean(x -> x.per, oscillatory_points)
        @info "Average period: $avg_period"

        avg_amplitude = mean(x -> x.amp, oscillatory_points)
        @info "Average amplitude: $avg_amplitude"


        # Store the results for the given fixed parameter combination
        results[(fixedpair[1].symbol, fixedpair[2].symbol)] = (volume = volume, avg_period = num_points, avg_amplitude, num_oscillatory_points = num_points)

        next!(loopprogress)
        display(loopprogress)
    end
    # Convert results to DataFrame
    # results_df = DataFrame(results)
    return results
end

# @trace define_parameter_constraints()
ic_constraints = define_initialcondition_constraints()
# @code_warntype GAProblem(param_constraints, fullprob)
ga_problem = GAProblem(ic_constraints, fullprob)
@profile run_GA(ga_problem; population_size=100) #? Vector of oscillatory points

Profile.print()
Profile.print(sortedby=:overhead, maxdepth=7)
Profile.clear()
@profile for i in 1:1000
    CostFunction(sol)
end
Profile.print()
@btime CostFunction($sol)
@code_warntype CostFunction(sol)




results = reachability_analysis(ic_constraints, fullprob)

eval_param_fitness(p,fullprob)

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

function testfunc(x)
    x+"hi"
    x[1] + x[2]
end

testfunc(1)


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