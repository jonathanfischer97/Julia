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
function make_ODEProb()#{Vector{Float64}}{Tuple{Float64, Float64}}
    #? Parameter list
    # psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
    #         :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
    #         :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    # p::Vector{Float64} = [x[2] for x in psym]
    p = [0.009433439939827041, 2.3550169939427845, 832.7213093872278, 12.993995997539924, 6.150972501791291, 
            1.3481451097940793, 0.006201726090609513, 0.006277294665474662, 0.9250191811994848, 57.36471615394549, 
            0.04411989797898752, 42.288085868394326, 3631.050539219606]    
    #? Initial condition list
    # usym = [:L => 0.0, :K => 0.5, :P => 0.3, :A => 2.0, :Lp => 3.0, :LpA => 0.0, :LK => 0.0, 
    #         :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
    #         :AKL => 0.0, :APLp => 0.0]
    # u0::Vector{Float64} = [x[2] for x in usym]
    u0 = [0.0, 0.5, 0.3, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    return ODEProblem(fullrn, u0, tspan, p) #TODO fix type stability
    # sol = solve(fullprob, saveat=0.1, save_idxs=1) #TODO fix type stability

    #? Plot the results
    # plot(sol)
end
@code_warntype make_ODEProb()
fullprob = make_ODEProb()

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

#TODO Testing type stability
fitness_function = make_fitness_function_with_fixed_inputs(eval_ic_fitness, fullprob, ic_constraints.ranges[1:2], (1, 2))
@code_warntype fitness_function(rand(14))
@code_warntype fitness_function(u0[1:end-2])

@benchmark solve($fullprob, saveat=0.1, save_idxs=1)
@code_warntype ODEProblem(fullrn, u0, tspan, p)
@code_warntype solve(fullprob, saveat=0.1, save_idxs=1)
@code_warntype eval_ic_fitness(rand(14), fullprob)

#<ChatGPT version revised
"""Monte carlo sampling to estimate bounding volume out of total volume"""
function monte_carlo_volume(points::Vector{Vector{Float64}}, var_constraints::Vector{ConstraintRange}, n_samples::Int = 10000)
    # Get the number of points and the dimension of each point
    n_points, dim = length(points), length(points[1])

    # Calculate the total volume using the parameter constraint ranges
    total_volume = prod([conrange.max - conrange.min for conrange in var_constraints])
    @info "Total solution volume: $total_volume"

    # Generate random points within the bounds specified by param_values
    random_points = generate_population(var_constraints, n_samples)

    # Calculate maximum distance threshold for points to be inside the hull
    max_dist_threshold = maximum([norm(points[i] .- points[j]) for i in 1:n_points for j in 1:n_points]) * 2

    # Initialize a counter for points inside the hull
    in_hull = zeros(Int64, Threads.nthreads())

    # For each random point, calculate the distance to each of the input points and check if it's in the hull
    Threads.@threads for i in 1:n_samples
        min_distance = minimum([norm(points[j] .- random_points[i]) for j in 1:n_points])

        # Check if the point is within the maximum distance threshold
        if min_distance <= max_dist_threshold
            in_hull[Threads.threadid()] += 1
        end
    end

    # Reduce the partial results into a final result
    total_in_hull = sum(in_hull)

    estimated_volume = total_volume * (total_in_hull / n_samples)

    return estimated_volume
end
#>End


generate_population(var_constraints,10000)
testvol = monte_carlo_volume([rand(2) for i in 1:10000], var_constraints)


using LazySets, Polyhedra, CDDLib


# Define the function to compute the result
function compute_result(points::Vector{Vector{Float64}}, var_constraints::ConstraintType)
    n = length(points)
    S = Vector{LazySet}(undef, n)

    @threads for i in 1:n
        P = Singleton(points[i])
        Q = HalfSpace(constraints[i].range[1].min, constraints[i].range[1].max)
        R = Hyperrectangle(constraints[i].range[2].min, constraints[i].range[2].max)
        S[i] = MinkowskiSum(Q, R)
    end

    # Compute the result
    result = Vector{Float64}(undef, n)
    @threads for i in 1:n
        result[i] = ρ(points[i], S[i])
    end

    return result
end
# Call the function
result = compute_result(points, constraints)


const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))


"""Compare the allowed solution space for when each pair of inputs is fixed"""
function reachability_analysis(constraints::ConstraintType, prob::ODEProblem) 

    #Get all combinations of fixed pairs
    fixed_pair_combos = combinations(constraints.ranges, 2)

    loopprogress = Progress(length(collect(fixed_pair_combos)), desc ="Looping thru fixed pairs: " , color=:blue)

    # Make a results dictionary where fixedpair => (volume, avg_period, avg_amplitude, num_oscillatory_points)
    results = Dict{Tuple{Symbol, Symbol}, NamedTuple{(:volume, :avg_period, :avg_amplitude, :num_oscillatory_points), Tuple{Float64, Float64, Float64, Int}}}() 

    # reachplot = barhist(title="Reachability Analysis", xlabel="Volume", ylabel="Average Period", zlabel="Average Amplitude", legend=:bottomright)
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

        # Run the optimization function to get the evaluated points
        oscillatory_points = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 8) #TODO: outputting the same number of points for multiple pairs
    
        # Compute convex hull volume of the oscillatory region in parameter space
        volume = monte_carlo_volume([x.ind for x in oscillatory_points], variable_constraints.ranges) #TODO: computing the same volume for multiple pairs
        @info "Volume of oscillatory region: $volume"
 

        # Calculate the number of oscillatory points
        num_points = length(oscillatory_points)
        @info "Number of oscillatory points: $num_points"

        # Calculate the average period and amplitude
        avg_period = mean(x -> x.per, oscillatory_points)
        @info "Average period: $avg_period"

        avg_amplitude = mean(x -> x.amp, oscillatory_points)
        @info "Average amplitude: $avg_amplitude"


        # Store the results for the given fixed parameter combination
        results[(fixedpair[1].symbol, fixedpair[2].symbol)] = (volume = volume, avg_period = avg_period, avg_amplitude, num_oscillatory_points = num_points)

        next!(loopprogress)
        # return results
    end
    # Convert results to DataFrame
    # results_df = DataFrame(results)
    return results
end

param_constraints = define_parameter_constraints()
param_reach_results = reachability_analysis(param_constraints, fullprob)

ic_constraints = define_initialcondition_constraints()
reach_results = reachability_analysis(ic_constraints, fullprob)




points = [x.ind for x in reach_results]
chull = convexhull(points...)
planarhull = Polyhedra.planar_hull(chull)
poly = polyhedron(planarhull)
poly = VPolygon(chull)
poly_area = area(poly)

# Make constraint hull
con_hull = convex_hull([[0.1,0.1],[0.1,10.0] ,[5.0,10.0],[5.0,0.1]])
constraint_poly = VPolygon(con_hull)
constraint_poly_area = area(constraint_poly)


plot(poly, ratio=:equal, legend=false, framestyle=:box, size=(500,500))
plot!([Singleton(p) for p in points[1:10:end]])

# Plot the constraint hull
plot!(constraint_poly, ratio=:equal, legend=false)



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