begin 
    using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using DiffEqCallbacks
    using Statistics
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter

    # using JET

    using LinearAlgebra

    using MultivariateStats
    using Clustering

    # using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

    # push!(LOAD_PATH, "/home/local/WIN/jfisch27/Desktop/Julia/UTILITIES")
    # using .OscTools

    using OscTools


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    FFTW.set_num_threads(18)
end






ogprobjac = make_ODE_problem();

callbackfunc = TerminateSteadyState(1e2)

sol = solve(ogprobjac, AutoTsit5(Rodas5P()), saveat=0.1, save_idxs=[6, 9, 10, 11, 12, 15, 16])
plot(sol)

cbsol = solve(ogprobjac, AutoTsit5(Rodas5P()), saveat=0.1, save_idxs=[6, 9, 10, 11, 12, 15, 16], callback=callbackfunc)
plot(cbsol)



paramconstraints = ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

#* calculate the minimum tspan in order to capture largest period oscillations
minimum_tspan = 1/(0.1 * 0.1) #multiply slowest forward rate by lowest concentration and take reciprocal to get seconds 

allconstraints = AllConstraints()
gaproblem = GAProblem(allconstraints, ogprobjac)
# @report_opt GAProblem(allconstraints, ogprobjac)

# @code_warntype GAProblem(allconstraints, ogprobjac)


initial_population = generate_population(allconstraints, 50000)
ga_results = run_GA(gaproblem, initial_population; iterations = 5)
ga_df = make_ga_dataframe(ga_results, allconstraints)
desc_df = describe(ga_df)
show(desc_df, allrows=true)
save_to_csv(ga_results, allconstraints, "adaptive_tspan_test.csv")

function test(gaproblem)
    Random.seed!(1234)
    initial_population = generate_empty_population(gaproblem.constraints, 3000)
    population = generate_population_stratified!(initial_population, gaproblem.constraints,10)
    run_GA(gaproblem, population; iterations = 5)
end

ga_results = test(gaproblem)

@btime test($gaproblem)

@bprofile test(gaproblem)

#* testing GA output for when fixed K=0.1, P=0.1, A=0.1, DF=100
ogprobjac = make_ODE_problem(100.);

gaprob = GAProblem(ode_problem = ogprobjac)
gaconstraints = gaprob.constraints

gaconstraints.DF.isfixed = true
gaconstraints.DF.fixed_value = 1000.0

gaconstraints.K.isfixed = true
gaconstraints.K.fixed_value = 0.01

gaconstraints.P.isfixed = true
gaconstraints.P.fixed_value = 0.1

gaconstraints.A.isfixed = true
gaconstraints.A.fixed_value = 0.1

gaprob
Random.seed!(1234)

initial_population = generate_population(gaconstraints, 10000)

ga_results = run_GA(gaprob, initial_population; iterations = 5)

#! population generation testing ##########
Random.seed!(1234)
loguniform_pop = generate_population(gaproblem.constraints, 5000)

initial_population = generate_empty_population(gaproblem.constraints, 5000)
stratified_pop = generate_population_stratified!(initial_population, gaproblem.constraints,3)

loguniform_df = make_pop_dataframe(loguniform_pop, gaproblem.constraints)
stratified_df = make_pop_dataframe(stratified_pop, gaproblem.constraints)


loguniform_stats = describe(loguniform_df)
show(loguniform_stats, allrows=true)
filter!(row -> row[1] != :DF, loguniform_stats)
CSV.write("loguniform_stats.csv", loguniform_df)

stratified_stats = describe(stratified_df)
show(stratified_stats, allrows=true)
filter!(row -> row[1] != :DF, stratified_stats)
CSV.write("stratified_stats.csv", stratified_df)

using StatsPlots

function plot_boxplots(df1::DataFrame, df2::DataFrame, title1::String, title2::String)
    p = plot(layout=(1,2), size=(1000, 500))
    @df df1 bar!(p[1], cols(16), title=title1, label = "P", yscale=:log10)
    @df df2 bar!(p[2], cols(16), title=title2, label = "P", yscale=:log10)
    display(p)
end


function plot_scatter(df1::DataFrame, df2::DataFrame)
    p = plot(size=(1000, 500), legend=:topright)
    for col in names(df1)
        if col != :variable
            scatter!(p, ["Min", "Max", "Median"], [df1[col, :min], df1[col, :max], df1[col, :median]], label="$col - Method 1")
            scatter!(p, ["Min", "Max", "Median"], [df2[col, :min], df2[col, :max], df2[col, :median]], label="$col - Method 2", shape=:diamond)
        end
    end
    xlabel!(p, "Statistic")
    ylabel!(p, "Value")
    title!(p, "Comparison of Min, Max, and Median Values")
    display(p)
end


plot_boxplots(loguniform_df, stratified_df, "Loguniform Population", "Stratified Population")



function shannon_diversity(population::Vector{Vector{Float64}}, nbins::Int)
    n_params = length(population[1])
    shannon_indices = zeros(n_params)
    
    for i in 1:n_params
        param_values = [indiv[i] for indiv in population]
        counts = fit(Histogram, param_values, nbins:nbins).weights
        p = counts ./ sum(counts)
        shannon_indices[i] = -sum(p .* log.(p))
    end
    
    return shannon_indices
end




