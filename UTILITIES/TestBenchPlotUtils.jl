#< Plotting utilities for testing
# """
# Plot the solution of an ODEProblem, `prob`, with respect to the variables `vars`.
# """
# function plotsol(prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), title = "")
#         sol = solve(prob, Rosenbrock23(), saveat=0.1, save_idxs=vars)
#         p = plot(sol, title = title, xlabel = "Time (s)", ylabel = "Concentration (µM)", lw = 2, size = (1000, 600))

#         display(p)
#         return p
# end


"""
Plot the solution from a row of the DataFrame
"""
function plotsol(sol::ODESolution; vars::Vector{Int} = collect(1:length(prob.u0)), fitidx = 1)

        cost, per, amp = CostFunction(sol)
        p = plot(sol, title = "Fit = $cost\nPeriod = $per s", xlabel = "Time (s)", ylabel = "Concentration (µM)", lw = 2, size = (1000, 600))
        # annotate!(p, (0, 0), text("Period = $per\nAmplitude = $amp", :left, 10))

        # display(p)
        return p    
end


"""
Plot the FFT of a solution from a row of the DataFrame
"""
function plotfft(sol::ODESolution; vars::Vector{Int} = collect(1:length(prob.u0)))

        normsol = normalize_time_series!(sol[1,:])
        # normsol = sol[1,:]
        solfft = getFrequencies(normsol)
        # fft_peakindexes, fft_peakvals = findmaxima(solfft,1) #* get the indexes of the peaks in the fft
        fft_peakindexes, peakprops = findpeaks1d(solfft; height = 1e-3, distance = 2) #* get the indexes of the peaks in the fft
        fft_peakvals = peakprops["peak_heights"]
        diffs = getDif(fft_peakvals)
        standevs = getSTD(fft_peakindexes, solfft)
        # cost, per, amp = CostFunction(sol)
        p1 = plot(solfft, title = "getDif: $(diffs)", xlabel = "Frequency (Hz)", ylabel = "Amplitude", lw = 2, xlims = (0, fft_peakindexes[end]), label="")
        peaklabels = [text("$(round.(val; digits=4))", :bottom, 8) for val in fft_peakvals]
        scatter!(p1, fft_peakindexes, fft_peakvals, text = peaklabels,
                         label = "Peakvals: $(round.(fft_peakvals; digits=4))\ngetSTD: $(standevs)\nSumFitness: $(-standevs-diffs)", color = :red, markersize = 5)
        p2 = plot(solfft, xlabel = "Frequency (Hz)", ylabel = "Amplitude", lw = 2, xlims = (0, 100), label="")
        scatter!(p2, fft_peakindexes, fft_peakvals, text = peaklabels, color = :red, markersize = 5)
        
        # display(p)
        return plot(p1, p2, layout = (2,1), size = (1000, 600))
end


"""
Plot both the solution and the FFT of a solution from a row of the DataFrame
"""
function plotboth(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), fitidx = 1)
        reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row]) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))])
        sol = solve(reprob, saveat=0.1, save_idxs=vars[1:5])

        solplot = plotsol(sol; vars)
        fftplot = plotfft(sol; vars)

        bothplot = plot(solplot, fftplot, layout = (2,1), size = (1000, 600))
        display(bothplot)
        return bothplot
end