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
function plotsol(sol::ODESolution; title = "", vars::Vector{Int} = collect(1:length(prob.u0)))

        # cost, per, amp = CostFunction(sol)
        p = plot(sol, idxs = vars, title = title, xlabel = "Time (s)", ylabel = "Concentration (µM)", lw = 2, size = (1000, 600))
        # annotate!(p, (0, 0), text("Period = $per\nAmplitude = $amp", :left, 10))

        # display(p)
        return p    
end


"""
Plot the FFT of a solution from a row of the DataFrame
"""
function plotfft(sol::ODESolution; fitidx=1)

        normsol = normalize_time_series!(sol[fitidx,:])
        # normsol = sol[1,:]
        solfft = getFrequencies(normsol)
        # fft_peakindexes, fft_peakvals = findmaxima(solfft,1) #* get the indexes of the peaks in the fft
        fft_peakindexes, peakprops = findpeaks1d(solfft; height = 1e-3, distance = 2) #* get the indexes of the peaks in the fft
        fft_peakvals = peakprops["peak_heights"]


        diffs = round(getDif(fft_peakvals); digits=4)
        standevs = round(getSTD(fft_peakindexes, solfft);digits=4)
        # cost, per, amp = CostFunction(sol)
        p1 = plot(solfft, title = "getDif: $(diffs)", xlabel = "Frequency (Hz)", ylabel = "Amplitude", lw = 2, xlims = (0, max(length(solfft),fft_peakindexes[end]+50)), label="", titlefontsize = 16, titlefontcolor = :green)
        peaklabels = [text("$(round.(val; digits=4))", :bottom, 10) for val in fft_peakvals]
        scatter!(p1, fft_peakindexes, fft_peakvals, text = peaklabels, label = "", color = :red, markersize = 5)

        window = 2
        maxpeak_idx = fft_peakindexes[argmax(fft_peakvals)]
        stdlines = [maxpeak_idx - window, maxpeak_idx + window]

        vline!(p1, stdlines, color = :blue, label = "")

        
        p2 = plot(solfft, title = "getSTD: $(standevs)", xlabel = "Frequency (Hz)", lw = 2, xlims = (0, 100), label="", titlefontsize = 16, titlefontcolor = :red)
        scatter!(p2, fft_peakindexes, fft_peakvals, text = peaklabels, color = :red, markersize = 5, label="")
        vline!(p2, stdlines, color = :blue, label = "")
        
        # display(p)
        return plot(p1, p2, layout = (1,2), size = (1000, 600))
end


"""
Plot both the solution and the FFT of a solution from a row of the DataFrame
"""
function plotboth(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), fitidx = 1)
        reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row]) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))])
        sol = solve(reprob, saveat=0.1, save_idxs=vars)
        cost, per, amp = CostFunction(sol; idx = fitidx)

        solplot = plotsol(sol; vars=vars[1:5])
        fftplot = plotfft(sol)

        bothplot = plot(solplot, fftplot, suptitle ="Fit: $(round(cost;digits=4))\nPeriod: $(round(per;digits=4)) s" , 
                        titlefont = :bold, layout = (2,1), size = (1000, 800))
        display(bothplot)
        return bothplot
end