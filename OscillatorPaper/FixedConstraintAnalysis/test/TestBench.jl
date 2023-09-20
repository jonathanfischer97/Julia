begin 
    using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using DiffEqCallbacks
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
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)



    using OscTools 


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    FFTW.set_num_threads(18)
end



#* Clustering test 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_10000/4FixedICRawSets/DF=1000.0")



