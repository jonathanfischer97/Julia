begin 
    # using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    # using DiffEqCallbacks
    using Statistics

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
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)



    using OscTools 


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW
    FFTW.set_num_threads(18)
end

# function test4fixedGA()
#     #* Set up the default GA problem
#     ga_problem = GAProblem()

#     #* Fixed some constraints 
#     set_fixed_constraints!(ga_problem.constraints, [:DF, :K, :P, :A])

#     #* Assign the fixed values 
#     set_fixed_values!(ga_problem.constraints, 1000., 1.0, 1.0, 3.16)

#     #* Set seed 
#     Random.seed!(1234)

#     #* Generate the initial population
#     population = generate_population(ga_problem.constraints, 10000)

#     #* Run the GA
#     run_GA(ga_problem, population)
# end

mutable struct EzraConstraints <: ConstraintSet
    ka1::ConstraintRange
    kb1::ConstraintRange
    ka4::ConstraintRange
    ka7::ConstraintRange
    DF::ConstraintRange

    L::ConstraintRange
    K::ConstraintRange
    P::ConstraintRange
    A::ConstraintRange
end

function EzraConstraints(; karange = (1e-2, 1.), kbrange = (1e-3, 1e2), dfrange = (1e4, 1e5),
                                Lrange = (1., 1e1), Krange = (1e-1, 1.), Prange = (1e-2, 1e-1), Arange = (1., 1e1))
    #* Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return EzraConstraints(
        ConstraintRange(name = :ka1, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb1, min = kb_min, max = kb_max),
        ConstraintRange(name = :ka4, min = ka_min, max = ka_max),
        ConstraintRange(name = :ka7, min = ka_min, max = ka_max),
        ConstraintRange(name = :DF, min = df_min, max = df_max),

        ConstraintRange(name = :L, min = Lrange[1], max = Lrange[2]),
        ConstraintRange(name = :K, min = Krange[1], max = Krange[2]),
        ConstraintRange(name = :P, min = Prange[1], max = Prange[2]),
        ConstraintRange(name = :A, min = Arange[1], max = Arange[2])
    )
end



OscTools.make_fitness_function_threaded(constraints::EzraConstraints, ode_problem::ODEProblem) = make_fitness_function_threaded(constraints, ode_problem, eval_ezra_fitness)


"""Evaluate the fitness of an individual with new initial conditions and new parameters"""
function eval_ezra_fitness(inputs::Vector{Float64}, prob::OP; idx::Vector{Int} = [6, 9, 10, 11, 12, 15, 16]) where {OP <: ODEProblem}
    newu = @view inputs[6:end]
    
    #* Allocate parameter array
    newp = zeros(13)
    newp[[1,2,8,10,13]] .= @view inputs[1:5] #assign sampled free parameters from input (ka1, kb1, ka4, ka7, df)
    newp[[4, 5, 6, 7, 12]] .= SA[0.7e-3, 2e-3, 0.118, 0.0609, 85.3] #assign known experimental parameter values (ka2, kb2, ka3, kb3, kcat7)
    
    ka1, kb1, ka4, ka7, df = inputs[1:5]

    #! Known experimental values to calculate other parameters
    Km1exp = 15.0
    Kd4exp =29.0
    Km7exp = 39.0
    kcat7exp = 85.3

    #! Calculate other parameters
    newp[3] = Km1exp * ka1 - kb1 #kcat1
    newp[9] = Kd4exp * ka4 #kb4
    newp[11] = Km7exp * ka7 - kcat7exp #kb7

    # @show newp
    # @show newu

    #! Check if any values in newp are negative, return [0.0, 0.0, 0.0] if so
    if any(newp .< 0.0)
        # @info "Negative parameter value, returning [0.0, 0.0, 0.0]"
        return [0.0, 0.0, 0.0]
    end
    

    #* remake with new initial conditions and new parameters
    new_prob = remake(prob; p = newp, u0= newu)
    return solve_for_fitness_peramp(new_prob, idx)
end


#* Make new GAProblem with Ezra's constraints
ezra_gaprob = GAProblem(constraints = EzraConstraints())

ez_pop = generate_population(ezra_gaprob.constraints, 300000)

ez_fitness_function = make_fitness_function_threaded(ezra_gaprob.constraints, make_ODE_problem(), eval_ezra_fitness)

ez_fitness_function(ez_pop[1])

eval_ezra_fitness(ez_pop[2], make_ODE_problem())

ezra_garesults = run_GA(ezra_gaprob, ez_pop)








#* Make figure for Maggie from Ezra's model 

"""
Converts row in oscillatory solution dataframe to an ODESolution
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
"""
function entryToSol(df, row)
    ka7Range = 10.0 .^((2:5) ./ 5)
    Km1exp = 15.0
    ka2exp = 0.7*10^-3
    kb2exp = 2.0 * 10^-3
    ka3exp = 0.118
    kb3exp = 0.0609
    Kd4exp =29.0
    Km7exp = 39.0
    kcat7exp = 85.3


    currow = df[row,:]
    u0 = zeros(16)
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    dfest = currow[:df]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
    # return solve(remake(prob, u0=u0, tspan=(0.0, tspan), p=p), Rodas4(), abstol=1e-8, reltol=1e-12, saveat=0.1, save_idxs=1, maxiters=200 * tspan, verbose=false, callback=TerminateSteadyState(1e-8, 1e-12))
    return u0, p 
end


# ezdf = CSV.read("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/EzraSolutions.csv", DataFrame)
ezdf = CSV.read("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/LocalData/AllExpOsc.csv", DataFrame)


oprob = make_ODE_problem()

#! Loop through all rows in ezdf and solve for AP2 oscillations
fitvals = zeros(nrow(ezdf), 3)
for i in 1:nrow(ezdf)
    newu, newp = entryToSol(ezdf, i)
    newinput = [newp; newu]
    fitvals[i,:] .= eval_all_fitness(newinput, oprob; idx = [1])
end

ezdf[!, :Cost] = fitvals[:,1]
ezdf[!, :Period] = fitvals[:,2]
ezdf[!, :Amplitude] = fitvals[:,3]

slow_newu, slow_newp = entryToSol(ezdf, 5)

###################!
testp = slow_newp[[1,2,8,10,13]]
testu = slow_newu[[1,2,3,4]]

testinput = [testp; testu]
append!(testinput, zeros(12))

eval_ezra_fitness(testinput, make_ODE_problem())
#####################!



slowprob = remake(oprob, u0=slow_newu, p=slow_newp, tspan=(0.,2000.))
slowsol = solve(slowprob, Rodas5P(), saveat=0.1)

fast_newu, fast_newp = entryToSol(ezdf, 3)
fastprob = remake(oprob, u0=fast_newu, p=fast_newp, tspan=(0.,2000.))
fastsol = solve(fastprob, Rodas5P(), saveat=0.1)

using Plots

function apply_default_settings(p)
    plot!(p, lw=4, size=(1000, 600), dpi=200,
          bottom_margin=12px, left_margin=16px, top_margin=10px, right_margin=8px)
    return p
end

function plotsol_maggie(sol::ODESolution; title = "")

    #* Sum up all A in solution 
    Asol = sol[4,:] .+ sol[13,:] .+ sol[14, :]

    #* Sum up all A on membrane
    Amem = OscTools.calculate_Amem(sol)
    
    p = plot(sol, idxs = [1,5,2,3], title = title, xlabel = "Time (s)", ylabel = "Concentration (µM)",
            color = [:blue :orange :purple :gold], label = ["PIP" "PIP2" "PIP5K" "Synaptojanin"], alpha = 0.5, lw=3)

    # p = plot(xlabel = "Time (s)", ylabel = "Concentration (µM)")

    plot!(p, sol.t, Asol, label="AP2 in solution", ls = :dash, alpha=1.0, color=:gray, lw=3)
    plot!(p, sol.t, Amem, label = "AP2 on membrane", ls = :dash, alpha=1.0, color=:black, lw=3)

    return p |> apply_default_settings
end


plotsol_maggie(fastsol; title= "Fast")

import CairoMakie as cm 

publication_theme() = cm.Theme(
    fontsize=16,
    font="CMU Serif",
    cm.Axis=(
        xlabelsize=20,xlabelpadding=-5,
        xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1,
        yticksize=10, xticksize=10,
        xlabel="x", ylabel="y"
        ),
    cm.Legend=(framecolor=(:black, 0.5), bgcolor=(:white, 0.5)),
    cm.Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5),
)

function maggie_plot(slowsol, fastsol)
    fig = cm.Figure(;resolution = (1000, 600))

    ax = cm.Axis(fig[1,1]; title= "Membrane binding of AP2",  xlabel = "Time (s)", ylabel = "Concentration (µM)",
    titlesize = 25, xlabelsize = 22, ylabelsize = 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1, ytickalign=1)

    cm.lines!(ax, slowsol.t, OscTools.calculate_Amem(slowsol), color=:blue, linewidth=3, label = "Period: 9.5 s")

    cm.lines!(ax, fastsol.t, OscTools.calculate_Amem(fastsol), color=:tomato, linewidth=3, label = "Period: 91 s")

    cm.axislegend(ax)

    fig
end

maggie_plot(slowsol, fastsol)


f, ax, l1 = cm.lines(slowsol.t, OscTools.calculate_Amem(slowsol);
                                                    figure = (; resolution = (1000, 600)),

                                                    axis = (; title= "Membrane binding of AP2", xlabel = "Time (s)", ylabel = "Concentration (µM)", xlabelsize = 20, ylabelsize = 20, ticklabelsize = 20, legendfontsize = 20),

                                                    color=:blue, linewidth=3, label = "Period: 9.5 s")
cm.lines!(ax, fastsol.t, OscTools.calculate_Amem(fastsol), color=:tomato, linewidth=3, label = "Period: 91 s")
cm.axislegend(ax)
f