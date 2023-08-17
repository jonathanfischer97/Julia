#! Solve model for arbitrary oscillatory parameters and initial conditions
function make_ODEProb(psym::Vector{Pair{Symbol, Float64}} = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
                                                        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
                                                        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606], 
                                        usym::Vector{Pair{Symbol, Float64}} = [:L => 3.0, :K => 0.5, :P => 0.3, :A => 2.0, :Lp => 0.0, :LpA => 0.0, :LK => 0.0, 
                                                :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]
                                        ;tspan = (0., 100.), rn = make_fullrn() 
                                                )
    #? Parameter list
    p= [x[2] for x in psym]
 
    #? Initial condition list
    u0 = [x[2] for x in usym]

    #? Create ODE problem
    return ODEProblem(rn, u0, tspan, p) #TODO fix type stability
end

function make_ODEProb(p::Vector{Float64} = [1.4353350626245021e-5, 0.20705634097197367, 461.9447701768986, 0.6642156268695386, 15.902417897116093, 7.971443130885463e-5, 0.7016715934086192, 
                                                1.5736952400326606e-5, 0.10411397662131976, 0.0007148707925785353, 0.05915700148314913, 0.8831745113213687, 1500.0], 
                                                u0::Vector{Float64} = [0.0, 301.15000000000003, 240.92000000000004, 602.3000000000001, 903.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                                                ;tspan = (0., 100.), rn = make_fullrn()
                                        )
        #? Create ODE problem
        return ODEProblem(rn, u0, tspan, p) #TODO fix type stability
end
#> TODO ^ Deprecate the functions above in favor of default values within the rn



#< Plotting utilities
"""
Plot the solution of an ODEProblem, `prob`, with respect to the variables `vars`.
"""
function plotsol(prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), title = "")
        sol = solve(prob, Rosenbrock23(), save_idxs=vars)
        p = plot(sol, title = title, xlabel = "Time (s)", ylabel = "Concentration (µM)", lw = 2, size = (1000, 600))

        display(p)
        return p
end


"""
Plot the solution from a row of the DataFrame
"""
function plotsol(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)))

        reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row]) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))])
        plotsol(reprob; vars, title = "Period = $(df.per[row]) s; Fitness = $(df.fit[row])")        
end


"""
Plot the FFT of a solution from a row of the DataFrame
"""
function plotfft(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)))

        reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row]) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))])
        sol = solve(reprob, Rosenbrock23(), saveat = 0.1, save_idxs=vars)
        normsol = normalize_time_series!(sol[1,:])
        # normsol = sol[1,:]
        solfft = getFrequencies(normsol)
        # fft_peakindexes, fft_peakvals = findmaxima(solfft,5) #* get the indexes of the peaks in the fft
        fft_peakindexes, peakprops = findpeaks1d(solfft; height = 1e-2) #* get the indexes of the peaks in the fft
        fft_peakvals = peakprops["peak_heights"]
        diffs = getDif(fft_peakvals)
        standevs = getSTD(fft_peakindexes, solfft)
        cost, per, amp = CostFunction(sol)
        p = plot(solfft, title = "Fit = $(cost)", xlabel = "Frequency (Hz)", ylabel = "Amplitude", lw = 2, size = (1000, 600), xlims = (0, 100), label="")
        scatter!(p, fft_peakindexes, fft_peakvals, label = "getDif: $(diffs)\nPeakvals: $(round.(fft_peakvals; digits=7))\ngetSTD: $(standevs)\nSumFitness: $(-standevs-diffs)\nPeriod: $per", color = :red, markersize = 5)
        
        display(p)
        return p
end


"""
Plot both the solution and the FFT of a solution from a row of the DataFrame
"""
function plotboth(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)))
        solplot = plotsol(row, df, prob; vars)
        fftplot = plotfft(row, df, prob; vars)

        bothplot = plot(solplot, fftplot, layout = (2,1), size = (1000, 600))
        display(bothplot)
        return bothplot
end