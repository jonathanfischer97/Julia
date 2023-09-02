begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    # using Peaks
    # using FindPeaks1D
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

tspan = (0., 2000.0)
fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], tspan, [])

de = modelingtoolkitize(ogprob)

ogprobjac = ODEProblem(de, [], tspan, jac=true)


# @code_warntype make_fitness_function(eval_param_fitness, ogprobjac; fitidx = 4)

# @btime solve($ogprob, saveat = 0.1, save_idxs = 1)
# @btime solve($ogprobjac, saveat = 0.1, save_idxs = 1)

# newprob = remake(ogprob, p = ogprob.p .* 1.0)
# @btime solve($newprob, Rosenbrock23(), saveat = 0.1, save_idxs = 1)

# newprobjac = remake(ogprobjac, p = ogprobjac.p .* 1.0)
# @btime solve($newprobjac, Rosenbrock23(), saveat = 0.1, save_idxs = 1)



# new_u0 = ogprob.u0 .* 10
# ogprob = remake(ogprob, u0 = new_u0)
# # @benchmark solve($ogprob, saveat = 0.1, save_idxs = 1)

# @benchmark solve($ogprob, saveat=0.1, save_idxs = 1)
# ogsol = solve(ogprob, saveat=0.1, save_idxs = 1)
# peakidxs, props = findpeaks1d(ogsol[1,:]; height = 0.1)
# plot(ogsol)
# testfunc(ogprob) = solve(ogprob, saveat=0.1, save_idxs = 1)

# using Cthulhu
# using ProfileView
# descend_code_warntype(testfunc, (ODEProblem,))
# @code_warntype solve(ogprob, saveat=0.1, save_idxs = 1)
# # plot(ogsol)
# @code_warntype solve_for_fitness_peramp(ogprob)

# @benchmark solve_for_fitness_peramp($ogprob)

# @code_warntype CostFunction(ogsol)
# @benchmark CostFunction($ogsol)
# @benchmark getPerAmp($ogsol)



#* Optimization of parameters to produce data for CSV
param_constraints = define_parameter_constraints(;karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = define_initialcondition_constraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

# gaproblem = GAProblem(param_constraints, ogprobjac)

# garesults = run_GA(gaproblem; population_size = 20000, iterations = 5)


function fixedDF_fitness_function_maker(evalfunc::Function, prob::ODEProblem, fixedDF::Float64)
    let evalfunc = evalfunc, prob = prob, fixedDF = fixedDF
        function fitness_function(input::Vector{Float64})
            newprob = remake(prob, p = [input; fixedDF])
            return evalfunc(input, newprob)
        end
    end
end


#* Function loops through 4D grid of different initial conditions, letting all parameters be freely optimized, and saves the results to a csv file
function fixed_quadruplet_ic_searcher(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength = 4, fixedDF=1000.)
    #* get the ranges of the initial conditions
    icranges = [logrange(constraints.min, constraints.max, rangelength) for constraints in icconstraints.ranges]

    icnames = [constraints.name for constraints in icconstraints.ranges]

    #* filter out DF because it will be fixed
    filter!(x -> x.name != "DF", paramconstraints.ranges)

    num_rows = rangelength^length(icnames)


    icvals1 = Vector{Float64}(undef, num_rows)
    icvals2 = Vector{Float64}(undef, num_rows)
    icvals3 = Vector{Float64}(undef, num_rows)
    icvals4 = Vector{Float64}(undef, num_rows)
    num_oscillatory_points_array = Vector{Int}(undef, num_rows)
    average_periods = Vector{Float64}(undef, num_rows)
    maximum_periods = Vector{Float64}(undef, num_rows)
    minimum_periods = Vector{Float64}(undef, num_rows)
    average_amplitudes = Vector{Float64}(undef, num_rows)
    maximum_amplitudes = Vector{Float64}(undef, num_rows)
    minimum_amplitudes = Vector{Float64}(undef, num_rows)


    i = 1

    #* make progress bar 
    loopprogress = Progress(num_rows, desc ="Looping thru fixed ICs: " , color=:red)

    mainrawpath = mkpath("./OscillatorPaper/FigureGenerationScripts/4FixedICRawSets")


    #* loop through each ic range and run the GA on each set of initial conditions after remaking the problem with them
    for icval1 in icranges[1]
        for icval2 in icranges[2]
            for icval3 in icranges[3]
                for icval4 in icranges[4]
                    icvals = [icval1, icval2, icval3, icval4]
                    @info icvals

                    #* remake the problem with the new initial conditions
                    newprob = remake(prob, u0 = [icvals; prob.u0[5:end]])
                    
                    #* make new GA problem with new initial conditions
                    ga_problem = GAProblem(paramconstraints, newprob)

                    #* set seed for reproducibility
                    Random.seed!(1234)

                    #* close fitness function maker 
                    fitness_function_maker(evalfunc, prob) = fixedDF_fitness_function_maker(evalfunc, prob, fixedDF)

                    #* run the GA on the new problem
                    oscillatory_points_results = run_GA(ga_problem, fitness_function_maker; population_size = 20000, iterations = 5)

                    #* get the number of oscillatory points
                    num_oscillatory_points = length(oscillatory_points_results.population)

                    #* if there are no oscillatory points, save the results to the results_df and continue
                    if iszero(num_oscillatory_points)
                        icvals1[i] = icval1
                        icvals2[i] = icval2
                        icvals3[i] = icval3
                        icvals4[i] = icval4
                        num_oscillatory_points_array[i] = 0
                        average_periods[i] = NaN
                        maximum_periods[i] = NaN
                        minimum_periods[i] = NaN
                        average_amplitudes[i] = NaN
                        maximum_amplitudes[i] = NaN
                        minimum_amplitudes[i] = NaN
                    else
                        average_periods[i]::Float64 = mean(oscillatory_points_results.periods)
                        maximum_periods[i]::Float64 = maximum(oscillatory_points_results.periods; init=0.0)
                        minimum_periods[i]::Float64 = minimum(oscillatory_points_results.periods; init=0.0)

                        average_amplitudes[i]::Float64 = mean(oscillatory_points_results.amplitudes)
                        maximum_amplitudes[i]::Float64 = maximum(oscillatory_points_results.amplitudes; init=0.0)
                        minimum_amplitudes[i]::Float64 = minimum(oscillatory_points_results.amplitudes; init=0.0)
                        
                        #* save the results to the results_df
                        # results_df[i, :] = (icval1, icval2, icval3, icval4, num_oscillatory_points, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude)
                        icvals1[i] = icval1
                        icvals2[i] = icval2
                        icvals3[i] = icval3
                        icvals4[i] = icval4
                        num_oscillatory_points_array[i] = num_oscillatory_points
                        
                        #* make dataframe from oscillatory_points_results
                        oscillatory_points_df = make_ga_dataframe(oscillatory_points_results, newprob, fixedDF)


                        innerrawpath = mkpath(mainrawpath*"/$(round(icval1; digits = 2))_$(round(icval2;digits = 2))_$(round(icval3; digits=2))_$(round(icval4; digits=2))")

                        CSV.write(innerrawpath*"/DF=$(round(fixedDF)).csv", oscillatory_points_df)
                    end
                    next!(loopprogress)
                    i += 1
                end
            end
        end
    end
    results_df = DataFrame(icnames[1] => icvals1, icnames[2] => icvals2, icnames[3] => icvals3, icnames[4] => icvals4,
                            "num_oscillatory_points" => num_oscillatory_points_array, 
                            "average_period" => average_periods, "maximum_period"=>maximum_periods, "minimum_period"=>minimum_periods,
                            "average_amplitude" => average_amplitudes, "maximum_amplitude"=>maximum_amplitudes, "minimum_amplitude"=>minimum_amplitudes)
                            
    CSV.write("./OscillatorPaper/FigureGenerationScripts/4FixedICs_DF=$(round(fixedDF)).csv", results_df)
    return results_df                
end

df = fixed_quadruplet_ic_searcher(param_constraints, ic_constraints, ogprobjac; rangelength=4, fixedDF=100.)

df = fixed_quadruplet_ic_searcher(param_constraints, ic_constraints, ogprobjac; rangelength=4, fixedDF=10000.)




function loop_4fixedICs_thru_DFvals(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength = 4, DFrange = [100.,1000.,10000.])
    for DF in DFrange
        df = fixed_quadruplet_ic_searcher(paramconstraints, icconstraints, prob; rangelength=rangelength, fixedDF=DF)
    end
end

loop_4fixedICs_thru_DFvals(param_constraints, ic_constraints, ogprobjac; rangelength=4, DFrange = [100.,1000.,10000.])


# scatter3d(df.Kinase, df.Phosphatase, df.AP2, color = df.num_oscillatory_points, zlims = (0, 100), colorbar_title = "Number of Oscillatory Points", title = "Number of Oscillatory Points vs. Initial Conditions", xlabel = "L", ylabel = "K", zlabel = "P", markersize = df.maximum_period, markerstrokewidth = 0, colormap = :viridis, camera = (30, 30),
#             xscale = :log10, yscale = :log10)  



df100 = CSV.read("./OscillatorPaper/FigureGenerationScripts/4FixedICs_DF=100.0.csv", DataFrame)
df1000 = CSV.read("./OscillatorPaper/FigureGenerationScripts/4FixedICs_DF=1000.0.csv", DataFrame)
df10000 = CSV.read("./OscillatorPaper/FigureGenerationScripts/4FixedICs_DF=10000.0.csv", DataFrame)

maximum(df10000.num_oscillatory_points)

testdf = CSV.read("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/4FixedICRawSets/100.0_0.22_0.22_1.0/DF=10000.0.csv", DataFrame)
plot_everything(testdf, ogprob; setnum=17, label="DF=100")


using GLMakie; GLMakie.activate!()

function scatters_in_3D(df; fig = Figure(resolution=(1600, 1200)))
    pnames = names(df)[1:4]
    xlog = log10.(df[:, pnames[1]])
    ylog = log10.(df[:, pnames[2]])
    zlog = log10.(df[:, pnames[3]])

    x = df[:, pnames[1]]
    y = df[:, pnames[2]]
    z = df[:, pnames[3]]

    # Identify non-NaN indices and values
    nonan_indices = findall(!isnan, df[:, :average_period])
    # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
    nonan_numpoints = df[:, :num_oscillatory_points][nonan_indices]

    nonan_periods = df[:, :average_period][nonan_indices]

    # Normalize sizes for non-NaN values
    sizes = fill(0.1, size(df, 1))
    # sizes[nonan_indices] = ((nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes))) ./ 2 
    sizes[nonan_indices] = ((nonan_numpoints .- minimum(nonan_numpoints)) ./ (maximum(nonan_numpoints) - minimum(nonan_numpoints))) ./ 2 


    # Normalize periods for non-NaN values
    norm_periods = fill(NaN, size(df, 1))

    norm_periods[nonan_indices] = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods)) 

    # Create the figure and axis
    
    ax = Axis3(fig[1:3,3]; aspect=:data, perspectiveness=0.5, title="Fixed ICs Oscillatory Regions", xlabel = pnames[1], ylabel = pnames[2], zlabel = pnames[3])

    # Scatter plot for non-NaN values
    hm = meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
                        diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true)

    # meshscatter!(ax, xlog, ylog; markersize=sizes, marker= Rect3f(Vec3f(0.,0.,0.1)), ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
    #                     diffuse = Vec3f(0.0), specular = Vec3f(0.0), shininess = 0, ambient = Vec3f(0.0))

    # Scatter plot for NaN values in gray
    # nan_indices = findall(isnan, df[:, :average_period])
    # meshscatter!(ax3, x[nan_indices], y[nan_indices], z[nan_indices]; markersize=sizes[nan_indices], color=:gray)

    # Colorbar and labels
    # Colorbar(fig[2, 1], hm, label="Period (s)", height=Relative(2.0))
    colgap!(fig.layout, 6)
    # xlabel!(ax3, "log10(kb3)")
    # ylabel!(ax3, "log10(kb4)")
    # zlabel!(ax3, "log10(DF)")

    # Display plot
    fig
end

fig = scatters_in_3D(df100)

fig2 = scatters_in_3D(df1000; fig=fig)




function plot4fixedIC(dfs; fig = Figure(resolution=(1600, 1200)))
    for (i,df) in enumerate(dfs)

        pnames = names(df)[2:4]
        xlog = log10.(df[:, pnames[1]])
        ylog = log10.(df[:, pnames[2]])
        zlog = log10.(df[:, pnames[3]])
        


        # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
        numpoints = df[:, :num_oscillatory_points]

        periods = df[:, :average_period]

        # Normalize sizes for non-NaN values
        sizes = fill(0.1, size(df, 1))
        # sizes[nonan_indices] = ((nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes))) ./ 2 
        sizes = ((numpoints .- minimum(numpoints)) ./ (maximum(numpoints) - minimum(numpoints))) ./ 2 


        # Normalize periods for non-NaN values
        norm_periods = fill(NaN, size(df, 1))

        norm_periods = (periods .- minimum(periods)) ./ (maximum(periods) - minimum(periods)) 

        # Create the figure and axis
        
        ax = Axis3(fig[1,i]; aspect=:data, perspectiveness=0.5, title="Fixed ICs Oscillatory Regions", xlabel = pnames[1], ylabel = pnames[2], zlabel = pnames[3])

        # Scatter plot for non-NaN values
        meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
                            diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true)
    end

    # meshscatter!(ax, xlog, ylog; markersize=sizes, marker= Rect3f(Vec3f(0.,0.,0.1)), ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
    #                     diffuse = Vec3f(0.0), specular = Vec3f(0.0), shininess = 0, ambient = Vec3f(0.0))


    # Colorbar and labels
    Colorbar(fig[2, 1], label="Period (s)", height=Relative(2.0))
    colgap!(fig.layout, 6)
    # xlabel!(ax3, "log10(kb3)")
    # ylabel!(ax3, "log10(kb4)")
    # zlabel!(ax3, "log10(DF)")

    # Display plot
    fig
end

fig = plot4fixedIC([df100, df1000])








#* compare 1000 pop to 100000 pop 
#* fix false positives, reduce std window, and height thresholding for find minima
#* for 3d scatter plots, show another plot where size is function of amplitude, and make amplitude the relative amplitude (Amem/Asol )
#* Maybe make plot for each value of L 
#* Increase grid sampling 