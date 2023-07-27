begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
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
    # using Combinatorics
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

    include("../../UTILITIES/ODEProbMaker.jl")

    include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end

fullrn = make_fullrn()

ratioprob = ODEProblem(fullrn, [], (0.,1000.0), [])

ratiosol = solve(ratioprob, Rosenbrock23(), save_idxs = 1)
plot(ratiosol)



function getFreqvsDF_threaded(prob, DFs=1000:100:100000)
    freqs = Vector{Float64}(undef,length(collect(DFs)))
    amps = Vector{Float64}(undef,length(collect(DFs)))
    Threads.@threads for i in eachindex(DFs)
        pcopy = deepcopy(prob.p)
        pcopy[13] = DFs[i]
        newprob = remake(prob, p = pcopy)
        sol = solve(newprob, save_idxs = 1, saveat = 0.01)
        fit, period, amplitude = CostFunction(sol)
        
        if fit == 0.0
            freqs[i], amps[i] = 0.0, 0.0
        else 
            freqs[i], amps[i] = 60/period, amplitude
        end
    end
    return DataFrames(collect(DFs), freqs, amps)
end



DFs, freqs, amps = getFreqvsDF_threaded(ratioprob)
DF_df = DataFrame(DF = DFs, freq = freqs, amp = amps)
CSV.write("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/DFvsFreq.csv", DF_df)




function plot_freqamp(varrange, freqs, amps, varlabel= "DF (nm)")
    freqplot = plot(varrange,freqs, xscale=:log10, legend=false, xlabel=varlabel, ylabel="Frequency (Hz)")
    ampplot = plot(varrange,amps, xscale=:log10, legend=false, xlabel=varlabel, ylabel="Amplitude (uM)")
    freqvsamp_plot = scatter(freqs, amps, legend=false, xlabel="Frequency", ylabel="Amplitude (uM)", ma=0.4)
    display(plot(freqplot, ampplot, freqvsamp_plot, plot_title = "$(ogprob.u0[1:4])", layout=(3,1), size=(1000,1000)))
end

plot_freqamp(DFs, freqs, amps)
# freqplot = plot(DFs,freqs, xscale=:log10, legend=false, xlabel="DF (nm)", ylabel="Frequency (Hz)")
# ampplot = plot(DFs,amps, xscale=:log10, legend=false, xlabel="DF (nm)", ylabel="Amplitude (uM)")
# freqvsamp_plot = scatter(freqs, amps, legend=false, xlabel="Frequency", ylabel="Amplitude (uM)", ma=0.4)
# plot(freqplot, ampplot, freqvsamp_plot, plot_title = "$(testnerdssprob.u0[1:4])", layout=(3,1), size=(1000,1000))


function getFreqvsConc(nerdssprob, concid1, concid2, concs=0.01:0.01:100.0)
    prob = deepcopy(nerdssprob)
    freqs = Vector{Float64}(undef,length(collect(concs)))
    amps = Vector{Float64}(undef,length(collect(concs)))
    for (i,conc) in enumerate(concs)
        prob.u0[concid1] = conc
        prob.u0[concid2] = conc
        nerdsol = solve(prob, save_idxs = 1, saveat = 0.01)
        fit, period, amplitude = CostFunction(nerdsol)
        
        if fit == 0.0
            freqs[i], amps[i] = 0.0, 0.0
        else 
            freqs[i], amps[i] = 60/period, amplitude
        end
    end
    return collect(concs), freqs, amps
end

concs, freqs, amps = getFreqvsConc(ratioprob, 2,3)
conc_DF = DataFrame(conc = concs, freq = freqs, amp = amps)
CSV.write("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ConcvsFreq.csv", conc_DF)
plot_freqamp(concs, freqs, amps, "PIP + AP2 (uM)")

testnerdssprob.p[13] = 1000000
testnerdsol = solve(testnerdssprob, save_idxs = 1)
getPerAmp(testnerdsol)
CostFunction(testnerdsol)
plot(testnerdsol)








#! run GA to get low copy number oscillations for NERDSS
ic_constraints = define_initialcondition_constraints(lipidrange = (0.1,1.5), kinaserange = (0.01,1.5), phosphataserange = (0.01,1.5), ap2range = (0.01,1.5))
param_constraints = define_parameter_constraints(karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e3, 2e4))

ic_gaproblem = GAProblem(ic_constraints,ogprob)
param_gaproblem = GAProblem(param_constraints,ogprob)

nerdss_record = run_GA(ic_gaproblem)


extract_solution(row) = extract_solution(row, nerdss_record, ogprob; vars = [1], tspan = (0., 100.))
plotsol(row) = plotsol(row, nerdss_record, ogprob; vars = [1])

plotsol(6)




function calc_concentration_ratio(dfrow::DataFrameRow)
    ratio = dfrow.ind[1]/sum(dfrow.ind[2:end-1])
    return round(ratio,digits = 3)
end

function calc_concentration_ratio(u0::Vector{Float64})
    ratio = u0[1]/sum(u0[2:end-1])
    return round(ratio,digits = 3)
end

minper = argmin(nerdss_record.per)
minpersol = extract_solution(minper)
minrow = nerdss_record[minper,:]
minfreq = 60/getPerAmp(minpersol)[1]
CostFunction(minpersol)
minper_concentration = calc_concentration_ratio(minrow)
minplot = plot(minpersol, label="Lipid:Enzyme ratio = $minper_concentration", xlabel = "Time (s)", ylabel = "PIP (μM)", 
            title = "Frequency = $(round(minfreq,digits = 3)) min⁻¹",legend_font_pointsize=14, legend=:topright)


maxper = argmax(nerdss_record.per)
maxpersol = extract_solution(maxper)
maxrow = nerdss_record[maxper,:]
maxfreq = 60/getPerAmp(maxpersol)[1]
CostFunction(maxpersol)
maxper_concentration = calc_concentration_ratio(maxrow)
maxplot = plot(maxpersol, label="Lipid:Enzyme ratio = $maxper_concentration",xlabel = "Time (s)", ylabel = "PIP (μM)",
                color = :blue, title = "Frequency = $(round(maxfreq,digits = 3)) min⁻¹",legend_font_pointsize=14)


minmaxplot = plot(minplot,maxplot, layout = (2,1),bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px, size = (1000,800),titlefontsize=18, guidefont = 14)
savefig(minmaxplot, "OscillatorPaper/FigureGenerationScripts/ProgressReportFigures/frequency_compare_plot5_LOWDF.png")







"""Interpolate intermediate initial conditions between min and max frequency solutions"""
function interpolate_points(start_point, end_point, num_points)
    # Ensure start_point and end_point are of same length
    if length(start_point) != length(end_point)
        error("Start point and end point must have the same number of dimensions")
    end

    # Create a range for each dimension
    ranges = [range(start_point[i], stop = end_point[i], length = num_points) for i in eachindex(start_point)]

    # Use the ranges to create an array of points
    points = [ [r[i] for r in ranges] for i in 1:num_points]

    return points
end

maxfreqind = minrow.ind 
minfreqind = maxrow.ind

#< INTERPOLATE BETWEEN MIN AND MAX FREQUENCY SOLUTIONS
intermediate_ics = interpolate_points(maxfreqind, minfreqind, 4)
# intermediate_ics = interpolate_points([1.25, 1.125, 0.856591382462083, 1.5], [1.440358036171343, 0.8507514631909991, 0.27735186553210145, 0.29092614167154873], 4)
nerdss_copyarrays = [concentration_to_copy_number.(ic, 0.4) for ic in intermediate_ics] #* Ordered from max frequency to min frequency 


function make_tunerplot(intermediate_ics)
    tunerplot = plot(layout = (4,1), legend=:none, plot_title = "Tuning the frequency of the oscillator", size = (1000,800), plot_titlefontsize=18, guidefont = 14)
    color_palette = palette([:red, :blue], length(intermediate_ics)) # Generate a color gradient
    
    for (i, ic) in enumerate(intermediate_ics)
        interprob = remake(ogprob, u0 = [ic; zeros(12)], tspan = (0.,20.))
        intersol = solve(interprob, save_idxs=1)
        interfreq = 60/getPerAmp(intersol)[1]
        plot!(tunerplot[i], intersol, title= "Frequency = $(round(interfreq,digits = 3)) min⁻¹", titlefontcolor = color_palette[i], titlelocation = :right, xlabel = "Time (s)", ylabel = "PIP (μM)", color = color_palette[i], label = "Lipid:Enzyme ratio = $(calc_concentration_ratio(ic))")
    end

    # Share x-axis across subplots
    # plot!(tunerplot, subplot = (4,1), link = :x)

    # display(tunerplot)
    return tunerplot
end

tunerplot = make_tunerplot(intermediate_ics)

savefig(tunerplot, "OscillatorPaper/FigureGenerationScripts/Figures/tunerplot1.png")



function extrapolate_points(start_point, end_point, num_points, extrapolate_factor=2.0)
    if length(start_point) != length(end_point)
        error("Start point and end point must have the same number of dimensions")
    end

    ranges = [range(start_point[i] - (end_point[i] - start_point[i]) * (extrapolate_factor - 1) / 2, 
                    stop = end_point[i] + (end_point[i] - start_point[i]) * (extrapolate_factor - 1) / 2, 
                    length = num_points) for i in eachindex(start_point)]

    points = [ [r[i] for r in ranges] for i in 1:num_points]

    return points
end

intermediate_ics = interpolate_points(maxfreqind, minfreqind, 100)
extrapolated_ics = extrapolate_points(maxfreqind, minfreqind, 100, 2.0)

#* Plot trend of concentrations vs frequency
function plot_ic_vs_frequency(intermediate_ics)
    freqs = []
    for (i, ic) in enumerate(intermediate_ics)
        interprob = remake(ogprob, u0 = [ic; zeros(12)], tspan = (0.,20.))
        intersol = solve(interprob, save_idxs=1)
        interfreq = 60/getPerAmp(intersol)[1]
        push!(freqs, interfreq)
    end
    conc_ratios = [calc_concentration_ratio(ic) for ic in intermediate_ics]
    plot(conc_ratios, freqs, xlabel = "Lipid:Enzyme Ratio", ylabel = "Frequency (min⁻¹)", legend = false, title = "Frequency vs Lipid:Enzyme ratio", titlefontsize=18, guidefont = 14)
end

icfreqplot = plot_ic_vs_frequency(extrapolated_ics)
savefig(icfreqplot, "OscillatorPaper/FigureGenerationScripts/Figures/icratio_freq_trendline.png")







# newpsym = [:ka1 => convert_to_macrorate(0.027267670545107203)#5.453534109021441e-05 #ka1, 1
#     :kb1 => 0.0643048008980449 #kb1, 2
#     :kcat1 => 286.6382253193995 #kcat1, 3
#     :ka2 => convert_to_macrorate(1.6605390671738467) #0.0033210781343476934 #ka2, 4
#     :kb2 => 0.39569337786534897 #kb2, 5
#     :ka3 => convert_to_macrorate(0.04115484270564225) #ka3, 6
#     :kb3 => 0.5393197910059361 #kb3, 7
#     :ka4 => convert_to_macrorate(0.05448532670350475) #ka4, 8
#     :kb4 => 0.2897657637531564 #kb4, 9
#     :ka7 => convert_to_macrorate(0.19013582091039463) #ka7, 10
#     :kb7 => 0.0028126177315505618 #kb7, 11
#     :kcat7 => 1.2733781341040291 #kcat7, 12
#     :DF => 5500.] #DF, 13
# newp = [x[2] for x in newpsym]


# newu0 = [4.9855314928351095, 0.08222720105366571, 0.669026631395984, 3.621476302891451]
newu0 = [0.9889272736535439, 0.754312949899292, 0.7*0.21447637196917263, 0.13832807513023648]
newp = [0.13631431320575543, 258.39273861543745, 922.897626015191, 
            1.7222377782178393, 49.044199345966945, 
                3.88687247761903, 0.5093036370614539, 
                0.6602571034403275, 1.3961637922644583, 
                10.627756359174143, 0.01, 33.6202081364354, 
                1.2*49096.34604920212]
newprob = remake(ogprob; u0 = [newu0; zeros(12)], p = newp, tspan = (0.,100.0))
newsol = solve(newprob, Rosenbrock23(), save_idxs = 1)
period, amplitude = getPerAmp(newsol)
freq = 60/period
@info "DF = $(newp[end]) nm, Frequency = $freq 1/min"
plot(newsol, title = "$(calc_concentration_ratio(newu0))", label = "Frequency = $(round(freq,digits = 3)) min⁻¹\nAmplitude = $(round(amplitude,digits = 3)) uM")

let
    DFs, freqs, amps = getFreqvsDF_threaded(newprob)
    plot_freqamp(DFs, freqs, amps)
end





function getFreqvsGrid(nerdssprob, concid1, concid2, concrange1=range(0.1, stop=100.0, length=1000), concrange2=range(10000, stop=100000.0, length=1000))
    prob = deepcopy(nerdssprob)
    freqs = Matrix{Float64}(undef,length(collect(concrange1)), length(collect(concrange2)))
    amps = Matrix{Float64}(undef,length(collect(concrange1)), length(collect(concrange2)))
    gridprogress = Progress(length(collect(concrange1))*length(collect(concrange2)), dt = 0.1, desc="Evaluating frequency grid... ", color=:red)
    Threads.@threads for i in eachindex(concrange1)
        for j in eachindex(concrange2)
            newu0 = copy(prob.u0)
            newu0[concid1] = concrange1[i]
            newp = copy(prob.p)
            newp[concid2] = concrange2[j]
            newprob = remake(prob; u0 = newu0, p = newp, tspan = (0.,1000.0))
            nerdsol = solve(newprob, save_idxs = 1, saveat = 0.1)
            fit, period, amplitude = CostFunction(nerdsol)
            
            if fit == 0.0
                freqs[j,i], amps[j,i] = 0.0, 0.0
                next!(gridprogress)
            else 
                freqs[j,i], amps[j,i] = 60/period, amplitude
                next!(gridprogress)
            end
        end
    end
    return collect(concrange1), collect(concrange2), freqs, amps
end


DFs, concs, freqs, amps = getFreqvsGrid(ogprob, 4, 13, range(0.1, stop=100.0, length=500), range(10000, stop=100000.0, length=500))

surface(DFs, concs, freqs, xlabel = "AP2 (μM)", ylabel = "V/A (nm)", zlabel = "Frequency min⁻¹", title = "Frequency vs V/A and AP2", color = :vik)
savefig("freqvsgrid.png")
contour(DFs, concs, freqs, xlabel = "AP2 (μM)", ylabel = "V/A (nm)", zlabel = "Frequency min⁻¹", title = "Frequency vs V/A and AP2", color = :vik)







randprob = remake(nerdssprob; u0 = vcat(nerdss_record[end,:].ind,zeros(12)), tspan = (0.,1000.0))
randsol = solve(randprob, Rosenbrock23(), saveat = 0.1)
plot(randsol)

function make_tuneplot(prob, idx1, idx2)
    per_array = []
    amp_array = []
    nerdsol = solve(prob, Rosenbrock23(), saveat = 0.1)
    per, amp = getPerAmp(nerdsol)
    pl = plot(nerdsol, idxs = (0,idx1,idx2), title = "Tuning the amplitude through V/A", label = "VA = $(prob.p[end]) nm\nPeriod = $(round(per)) s\nAmplitude = $(round(amp)) uM\n")

    pcopy = copy(prob.p)
    for i in 1:2
        pcopy[end] = pcopy[end] * 1.5;
        @info "V/A = $(pcopy[end]) nm"
        new_nerdsol = solve(remake(prob; p = pcopy), Rosenbrock23(), saveat = 0.1);
            per, amp = getPerAmp(new_nerdsol)
            println("Period = $per")
            println("Amplitude = $amp")
            push!(per_array, per)
            push!(amp_array, amp)
        plot!(pl, new_nerdsol, idxs = (0,idx1,idx2),label = "V/A = $(pcopy[end]) nm\nPeriod = $(round(per)) s\nAmplitude = $(round(amp)) uM\n", ls = :dash)
    end
    display(pl)
    return pl
end

p1 = make_tuneplot(nerdssprob, 1, 10)
savefig(p1, "/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ProgressReportFigures/tuneplot1.png")

#* changing initial conditions to tune period
begin 
    per_array = []
    amp_array = []

    pl = plot(nerdsol, idxs = (0,1,3), xlabel = "PIP2 (uM)", ylabel = "Synaptojanin (µM)", title = "Tuning the period through V/A", label = "Period = $(getPerAmp(nerdsol)[1]) s")

    for i in 1:3
        u0[4] = u0[4] * 0.8;
        @info "V/A = $(u0[4]) nm"
        new_nerdsol = solve(remake(nerdssprob; u0 = u0), Rosenbrock23(), saveat = 0.1);
            per, amp = getPerAmp(new_nerdsol)
            println("Period = $per")
            println("Amplitude = $amp")
            push!(per_array, per)
            push!(amp_array, amp)
        plot!(pl, new_nerdsol, idxs = (0,1,3),label = "Period = $per s", ls = :dash)
    end
    u0 = [x[2] for x in usym]
    display(pl)
end
savefig(pl, "OscillatorPaper/FigureGenerationScripts/ProgressReportFigures/NERDSS_Tuning.png")

plot(nerdsol, xlabel = "PIP2 (uM)", ylabel = "Synaptojanin (µM)", title = "Tuning the period through V/A", label = "VA = $(p[end]) s")
getPerAmp(nerdsol)

findmaxima(nerdsol[1,:],1000)