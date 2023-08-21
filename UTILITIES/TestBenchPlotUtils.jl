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
Calculates the frequency per minute of the FFT vector in-place 
"""
function frequencies_per_minute!(t::Vector{Float64}, freq_values::Vector{Float64})
        # Calculate the time step between samples
        dt = t[2] - t[1]
        
        # Calculate the length of the original solution array
        N = 2 * length(freq_values)
        
        # Calculate the frequency step in Hz (per second)
        freq_step_per_second = 1 / (N * dt)
        
        # Convert the frequency step to per minute
        freq_step_per_minute = freq_step_per_second * 60
        
        # Update the frequency values in-place to frequencies per minute
        freq_values .= freq_values .* freq_step_per_minute
    end
    
    
    


"""
Plot the FFT of a solution from a row of the DataFrame
"""
function plotfft(sol::ODESolution; fitidx::Int=4)
        #* Trim first 10% of the solution array to avoid initial spikes
        tstart = cld(length(sol.t),10) 
        trimsol = sol[tstart:end] 

        # normsol = normalize_time_series!(trimsol[fitidx,:])
        # normsol = sol[1,:]
        solfft = getFrequencies(trimsol[fitidx,:])

        #* Get the frequencies per minute for x axis
        frequencies_per_minute!(sol.t, solfft)

        #* Normalize the FFT to have mean 0 and amplitude 1
        normalize_time_series!(solfft)
        # fft_peakindexes, fft_peakvals = findmaxima(solfft,1) #* get the indexes of the peaks in the fft
        fft_peakindexes, peakprops = findpeaks1d(solfft; height = 1e-3, distance = 2) #* get the indexes of the peaks in the fft
        
        #* If there are no peaks, return a plot with no peaks
        if isempty(fft_peakindexes)
                p1 = plot(solfft, title = "getDif: 0.0", xlabel = "Frequency (min⁻¹)", ylabel = "Amplitude", lw = 2, 
                                xlims = (0, 100), label="", titlefontsize = 18, titlefontcolor = :green)
                return p1
        else
                fft_peakvals = peakprops["peak_heights"]
                diffs = round(getDif(fft_peakvals); digits=4)
                standevs = round(getSTD(fft_peakindexes, solfft; window = 5);digits=4)


                p1 = plot(solfft, title = "getDif: $(diffs)", xlabel = "Frequency (min⁻¹)", ylabel = "Amplitude", lw = 2, 
                                xlims = (0, min(length(solfft),fft_peakindexes[end]+50)), ylims=(0.0,min(1.0, maximum(fft_peakvals)+0.25)), label="", titlefontsize = 18, titlefontcolor = :green)
                peaklabels = [text("$(round.(val; digits=4))", :bottom, 10) for val in fft_peakvals]
                scatter!(p1, fft_peakindexes, fft_peakvals, text = peaklabels, label = "", color = :red, markersize = 5)

                window = 1
                maxpeak_idx = fft_peakindexes[argmax(fft_peakvals)]
                stdlines = [maxpeak_idx - window, maxpeak_idx + window]

                
                p2 = plot(solfft, title = "getSTD: $(standevs)", xlabel = "Frequency (min⁻¹)", lw = 2, xlims = (max(0,maxpeak_idx-50), min(length(solfft),maxpeak_idx+50)), 
                                                ylims=(0.0,min(1.0, maximum(fft_peakvals)+0.25)),label="", titlefontsize = 18, titlefontcolor = :red)
                scatter!(p2, fft_peakindexes, fft_peakvals, text = peaklabels, color = :red, markersize = 5, label="")
                vline!(p2, stdlines, color = :blue, label = "")
                
                return plot(p1, p2, layout = (1,2), size = (1000, 600))
        end
end


s"""
Plot both the solution and the FFT of a solution from a row of the DataFrame
"""
function plotboth(row, df::DataFrame, prob::ODEProblem; vars::Vector{Int} = collect(1:length(prob.u0)), fitidx::Int = 4)

        reprob = length(df.ind[row]) > 4 ? remake(prob, p = df.ind[row]) : remake(prob, u0 = [df.ind[row]; zeros(length(prob.u0) - length(df.ind[row]))])
        sol = solve(reprob, Rosenbrock23(), saveat=0.1, save_idxs=vars)
        cost, per, amp = CostFunction(sol; idx = fitidx)

        solplot = plotsol(sol; vars=vars[1:5])
        fftplot = plotfft(sol; fitidx = fitidx)

        bothplot = plot(solplot, fftplot, plot_title ="Fit: $(round(cost;digits=4))\nPeriod: $(round(per;digits=4)) s" , 
                        plot_titlefontsize = 20, layout = (2,1), size = (1000, 800))
        display(bothplot)
        return bothplot
end

"""
Plot the solution and FFT of every row in the DataFrame
"""
function plot_everything(df::DataFrame, prob::ODEProblem; setnum::Int = 1, label = "", jump=10)
        fftw_threads = FFTW.get_num_threads()
        FFTW.set_num_threads(18)
        progbar = Progress(cld(nrow(df),jump); desc = "Plotting:")
        path = mkpath("OscillatorPaper/FigureGenerationScripts/TestbenchPlots/Set$(setnum)-$(label)")
        CSV.write(path*"/Set$(setnum)-$(label).csv", df)
    
        for i in 1:jump:nrow(df)
            p = plotboth(i, df, prob)
            savefig(p, path*"/plot_$(i).png")
            next!(progbar)
        end
        FFTW.set_num_threads(fftw_threads)
    end