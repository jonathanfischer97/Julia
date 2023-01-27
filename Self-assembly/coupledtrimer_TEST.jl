using Random 
using Plots 
using Catalyst
using DifferentialEquations
# using Combinatorics
using Colors
using ColorSchemes
# using ColorTypes
# using LinearAlgebra
using SpecialFunctions
using ProgressMeter
using Plots.PlotMeasures
# using Latexify

trimer_rn = @reaction_network trirn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*y,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 
    (ka8,kb8), LpA + T <--> X

    (ka1m, kb1m), X + B <--> XB 
    (ka1m, kb1m), X + C <--> XC 
    (ka1m, kb1m), X + D <--> XD 
    (ka4m, kb4m), XB + C <--> XBC
    (ka4m*y, kb4m), XB + XC <--> XXBC
    (ka4m, kb4m), XC + B <--> XBC
    (ka5m, kb5m), XB + D <--> XBD
    (ka5m*y, kb5m), XB + XD <--> XXBD
    (ka5m, kb5m), XD + B <--> XBD
    (ka6m, kb6m), XC + D <--> XCD
    (ka6m*y, kb6m), XC + XD <--> XXCD
    (ka6m, kb6m), XD + C <--> XCD 
    (ka7m, kb7m), XBC + D <--> XBCD 
    (ka7m, kb7m), XXBC + D <--> XXBCD
    (ka7m*y, kb7m), XXBC + XD <--> XXXBCD
    (ka8m, kb8m), XBD + C <--> XBCD 
    (ka8m, kb8m), XXBD + C <--> XXBCD 
    (ka8m*y, kb8m), XXBD + XC <--> XXXBCD
    (ka9m, kb9m), XCD + B <--> XBCD 
    (ka9m, kb9m), XXCD + B <--> XXBCD 
    (ka9m*y, kb9m), XXCD + XB <--> XXXBCD 
    (ka1m*y, kb1m), XBC + X <--> XXBC 
    (ka1m*y, kb1m), XBD + X <--> XXBD 
    (ka1m*y, kb1m), XCD + X <--> XXCD 
    (ka1m*y, kb1m), XBCD + X <--> XXBCD 
    (ka1m*y, kb1m), XXBCD + X <--> XXXBCD 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 ka8 kb8 ka1m kb1m ka4m kb4m ka5m kb5m ka6m kb6m ka7m kb7m ka8m kb8m ka9m kb9m y 

#parameter list
p = [14.766654049574964, 247.74197239200362, 500.0, 31.051821195245743, 158.0304365645891, 52.425931974326005, 500.0, 29.68965224106596, 
3.1517004441718446, 90.20783657199675, 364.05856621148223, 96.68089672060627, 100.0, 178.4206713875396, 63.913390412700274, 16.670218374002296, 
31.00096833128655, 187.53746678023433, 24.208313355774997, 299.70860654065797, 62.268232479997536, 500.0, 13.217117644835021, 69.69105864551203, 
79.85635210382746, 247.8070093520743, 6.034069612733998, 28.445138804106804, 0.7817627051993913/0.001]

pmap = (:ka1 => 14.766654049574964, :kb1 => 247.74197239200362, :kcat1 => 500.0, :ka2 => 31.051821195245743, :kb2 => 158.0304365645891, :ka3 => 52.425931974326005, :kb3 => 500.0, :ka4 => 29.68965224106596, :kb4 => 3.1517004441718446, :ka7 => 90.20783657199675, :kb7 => 364.05856621148223, :kcat7 => 96.68089672060627, :ka8 => 100.0, :kb8 => 178.4206713875396, :ka1m => 63.913390412700274, :kb1m => 16.670218374002296, :ka4m => 31.00096833128655, :kb4m => 187.53746678023433, :ka5m => 24.208313355774997, :kb5m => 299.70860654065797, :ka6m => 62.268232479997536, :kb6m => 500.0, :ka7m => 13.217117644835021, :kb7m => 69.69105864551203, :ka8m => 79.85635210382746, :kb8m => 247.8070093520743, :ka9m => 6.034069612733998, :kb9m => 28.445138804106804, :y => 0.7817627051993913/0.001)
#initial condition list
u0 = [0. , 0.2, 0. , 3. , 0.9, 0. , 0. , 0. , 0.3, 0. , 0. , 0. , 0.3,
1. , 0.5, 0., 0.5, 0.0 , 0.5 , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. ]


# pmap = [x => y for (x,y) in zip(parameters(trimer_rn),p)]
# pmap = symmap_to_varmap(trimer_rn, pmap)

umap = [x => y for (x,y) in zip(states(trimer_rn),u0)]
umap = symmap_to_varmap(trimer_rn, umap)

tspan = (0.,20.)

odesys = convert(ODESystem, trimer_rn; combinatoric_ratelaws=true)
oprob = ODEProblem(trimer_rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5(), saveat = 0.001)


## Time Series gif ##
function plotfig1!(plt, osol;tend=length(osol))

    tend = Int(round(tend))

    # First Plot
    plot!(plt,osol.t[1:tend], osol[1,1:tend], color = :blue, label = "Non-sticky lipid")
    plot!(plt, osol.t[1:tend], osol[4,1:tend], color = :orange, label = "Sticky lipid")
    plot!(plt, osol.t[1:tend], osol[5,1:tend], color = :red, label = "Adaptor")
    # plot!(plt, osol.t[1:tend], osol[6,1:tend], color = :green, label = "PIP2-AP2")
    plot!(plt, osol.t[1:tend], osol[7,1:tend], color = :gold, label = "Membrane-bound kinase")
    plot!(plt, osol.t[1:tend], osol[11,1:tend], color = :purple, label = "Membrane-bound phosphatase")
end



default(linewidth = 3.0)

function make_lipid_timeseries_gif(osol)
    tframe = range(1,stop=length(osol.t),length=100)
    gif1 = @animate for i in tframe
        # idx = Int(round(i))
        plt1 = plot(dpi = 200,size = (1000, 400), legend = :outertopright, bottom_margin = 12px, left_margin = 16px, framestyle=:box, legendfontsize=10, linewidth = 3.0, xlims = (0,20), ylims = (0,3.5))
        # ylims!(plt[1], (0, 3.5))
        xlabel!("Time (s)")
        ylabel!("Concentration (uM)")
        # title!("Trimerization Reaction")
        plotfig1!(plt1,osol;tend=i)
    end 
    gif(gif1, "lipid_timeseries.gif")
end


make_lipid_timeseries_gif(osol)


## SSA ## 
function copynumber(conc, volume) #converts concentration to copy number
    newvolume = volume/1e15 #converts liters to um^3
    newconc = conc/1e6 #converts umol to mol
    moles = newconc * newvolume #volume must be passed in um^-3
    moles * 6.023e23
end

function convert_rate(rate, volume) #converts rate from us^-1 to copy number per second
    rate1 = rate/0.602214076 ##conversion ratio from page 10 of NERDSS manual
    rate2 = rate1*1e6 ##convert from us to s
    newvolume = volume*1e9 #convert volume to nm3
    rate2/newvolume ##calculate copy numbers per second
end

function convert_on_rates(pmap)
    p = []
    for (k,v) in pmap
        if occursin("ka",string(k))
            push!(p,convert_rate(v, 0.5))
        else
            push!(p,v)
        end
    end
    return p
end



# next we create a discrete problem to encode that our species are integer valued:
dprob = DiscreteProblem(trimer_rn, copynumber.(u0, 0.5), tspan, convert_on_rates(pmap))

# now, we create a JumpProblem, and specify Gillespie's Direct Method as the solver:
jprob = JumpProblem(trimer_rn, dprob, Direct(), save_positions=(false,false))

# now, let's solve and plot the jump process:
sol = solve(jprob, SSAStepper(), saveat=0.001)



function make_lipid_timeseries_ssa_gif(osol)
    tframe = range(1,stop=length(osol.t),length=100)
    gif1 = @animate for i in tframe
        # idx = Int(round(i))
        plt1 = plot(dpi = 200,size = (1000, 400), legend = :outertopright, bottom_margin = 12px, left_margin = 16px, framestyle=:box, legendfontsize=10, linewidth = 3.0, xlims = (0,20), ylims = (0,1100))
        # ylims!(plt[1], (0, 3.5))
        xlabel!("Time (s)")
        ylabel!("Copy Number")
        # title!("Trimerization Reaction")
        plotfig1!(plt1,osol;tend=i)
    end 
    gif(gif1, "lipid_timeseries_ssa.gif")
end


function make_lipid_timeseries_subplots_gif(osol, sol)
    tframe = range(1,stop=length(osol.t),length=100)
    gif1 = @animate for i in tframe
        fig = plot(layout = grid(2,1), dpi = 200,size = (1000, 600), legend = false, bottom_margin = 12px, left_margin = 16px, framestyle=:box, legendfontsize=10, linewidth = 3.0, xlims = (0,20))
        # ylims!(plt[1], (0, 3.5))
        xlabel!(fig[2],"Time (s)")
        ylabel!(fig[1],"Concentration (uM)")
        ylabel!(fig[2],"Copy Number")
        ylims!(fig[1], (0, 3.5))
        ylims!(fig[2], (0, 1100))

        plotfig1!(fig[1],osol;tend=i)
        plotfig1!(fig[2],sol;tend=i)
    end 
    gif(gif1, "lipid_timeseries_subplots.gif")
end

make_lipid_timeseries_subplots_gif(osol, sol)







default(linewidth = 3.0, dpi = 200)


function find_idx(pmap, p)
    for (i, (k,v)) in enumerate(pmap)
        if string(k) == p
            return i
        end
    end
end

function make_2dphaseplot()
    p1 = plot(osol[2,:],osol[end,:], label = "V/A")

    pcopy = copy(p)
    idx = find_idx(pmap, "y")

    color_gradient = cgrad(:plasma, 5, categorical = true)

    for (i,c) in enumerate(color_gradient)
        pcopy[idx] = pcopy[idx]*1.1
        println(pcopy[idx])
        newpmap = [x => y for (x,y) in zip(parameters(trimer_rn),pcopy)]
        newpmap = symmap_to_varmap(trimer_rn, newpmap)
        newprob = ODEProblem(trimer_rn, umap, tspan, newpmap)
        newsol = solve(newprob, Tsit5())
        plot!(p1,newsol[2,:],newsol[end,:],c = c, label = "V/A * $(round(1.1^i;digits=1))")
    end

    scatter!((osol[2,1],osol[end,1]), label = "Start", markersize = 5, color = :green)
    scatter!((osol[2,end],osol[end,end]), label = "End", markersize = 5, color = :red)
    xlabel!("PIP5K Kinase (uM)")
    ylabel!("XXXBCD Trimer (uM)")
    display(p1)
end

make_2dphaseplot()




function make_2dphase_gif(osol, x, y, param, scalar) #make gif of 2D phase plot 
    tframe = range(2,stop=length(osol.t),length=500)
    tframevec = [Int(round(x)) for x in tframe]
    color_gradient = cgrad(:plasma, 6, categorical = true)
    pcopy = copy(p)
    idx = find_idx(pmap, param)
    xidx, yidx = find_idx(umap, x), find_idx(umap, y)

    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        pcopy[idx] = pcopy[idx]*scalar
        newprob = ODEProblem(trimer_rn, umap, tspan, pcopy)
        newsol = solve(newprob, Tsit5(), saveat = 0.001)
        newsolarray[i] = newsol
    end

    prog = Progress(length(tframevec), 1, "Making gif: ")
    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot(osol[xidx,1:f], osol[yidx,1:f], c = color_gradient[1], seriesalpha = alpharange, xlabel = x, ylabel = y, label = "", xlims = (0,0.25), ylims = (0,3.5))
        scatter!((osol[xidx,f],osol[yidx,f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = "")
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[xidx,1:f],nsol[yidx,1:f],c = c, seriesalpha = alpharange, label = "")
            scatter!((nsol[xidx,f],nsol[yidx,f]), markershape = :circle, markersize = 7, c=c, label = "")#param * " * $(round(1.1^i;digits=1))")
        end
        next!(prog)
    end  
    gif(anim, x*"-"*y*"_"*param*"_"*"2dphase.gif")
end

make_2dphase_gif(osol, "XXXBCD(t)", "Lp(t)", "y", 1.1)




function make_2dphase_gif_INITCONDITIONS(osol, x, y, species, scalar) #make gif of 2D phase plot 
    tframe = range(2,stop=length(osol.t),length=300)
    tframevec = [Int(round(x)) for x in tframe]
    color_gradient = cgrad(:plasma, 6, categorical = true)
    ucopy = copy(u0)
    idx = find_idx(umap, species)
    xidx, yidx = find_idx(umap, x), find_idx(umap, y)

    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        ucopy[idx] = ucopy[idx]*scalar
        newprob = ODEProblem(trimer_rn, ucopy, tspan, p)
        newsol = solve(newprob, Tsit5(), saveat = 0.001)
        newsolarray[i] = newsol
    end

    prog = Progress(length(tframevec), 1, "Making gif: ")
    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot(osol[xidx,1:f], osol[yidx,1:f], c = color_gradient[1], seriesalpha = alpharange, xlabel = x, ylabel = y, label = "", xlims = (0,0.25), ylims = (0,0.007))
        scatter!((osol[xidx,f],osol[yidx,f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = species)
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[xidx,1:f],nsol[yidx,1:f],c = c, seriesalpha = alpharange, label = "")
            scatter!((nsol[xidx,f],nsol[yidx,f]), markershape = :circle, markersize = 7, c=c, label = species * " * $(round(1.1^i;digits=1))")
        end
        next!(prog)
    end  
    gif(anim, x*"-"*y*"_"*species*"_"*"2dphaseINIT.gif")
end

make_2dphase_gif_INITCONDITIONS(osol, "XXXBCD(t)", "A(t)", "D(t)", 1.2)



function make_2dphase_gif_SSA(osol, x, y, param, scalar) #make gif of 2D phase plot 
    tframe = range(2,stop=length(osol.t),length=300)
    tframevec = [Int(round(x)) for x in tframe]
    color_gradient = cgrad(:plasma, 6, categorical = true)
    pcopy = copy(p)
    idx = find_idx(pmap, param)
    xidx, yidx = find_idx(umap, x), find_idx(umap, y)

    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        pcopy[idx] = pcopy[idx]*scalar
        newprob = ODEProblem(trimer_rn, umap, tspan, pcopy)
        newsol = solve(newprob, Tsit5(), saveat = 0.001)
        newsolarray[i] = newsol
    end

    prog = Progress(length(tframevec), 1, "Making gif: ")
    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot(osol[xidx,1:f], osol[yidx,1:f], c = color_gradient[1], seriesalpha = alpharange, xlabel = x, ylabel = y, label = "")
        scatter!((osol[xidx,f],osol[yidx,f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = param)
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[xidx,1:f],nsol[yidx,1:f],c = c, seriesalpha = alpharange, label = "")
            scatter!((nsol[xidx,f],nsol[yidx,f]), markershape = :circle, markersize = 7, c=c, label = param * " * $(round(1.1^i;digits=1))")
        end
        next!(prog)
    end  
    gif(anim, x*"-"*y*"_"*param*"_"*"2dphase.gif")
end

make_2dphase_gif(osol, "XXXBCD(t)", "A(t)", "ka2", 1.2)


















function make_3dphase_gif(osol) #make gif of #D phase plot 
    color_gradient = cgrad(:plasma, 6, categorical = true)
    pcopy = copy(p)
    idx = find_idx(pmap, "y")

    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        pcopy[idx] = pcopy[idx]*1.1
        newpmap = [x => y for (x,y) in zip(parameters(trimer_rn),pcopy)]
        newpmap = symmap_to_varmap(trimer_rn, newpmap)
        newprob = ODEProblem(trimer_rn, umap, tspan, newpmap)
        newsol = solve(newprob, Tsit5(), saveat = 0.01)
        newsolarray[i] = newsol
    end

    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot3d(osol[2,1:f], osol[end,1:f], osol[4,1:f], arrow = true, c = color_gradient[1], seriesalpha = alpharange, xlims=(0.12,0.2),ylims=(0,0.22), label = param, camera=(360k/100, 30))
        # scatter3d!((osol[2,f],osol[end,f],osol[4,1:f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = "V/A")
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[2,1:f],nsol[end,1:f],nsol[4,1:f], arrow = true, c = c, seriesalpha = alpharange, label = param * " * $(round(1.1^i;digits=1))")
            # scatter3d!((nsol[2,f],nsol[end,f],nsol[4,1:f]), markershape = :circle, markersize = 7, c=c, label = "V/A * $(round(1.1^i;digits=1))")

        end
    end  
    gif(anim, "3dphase.gif")
end

make_3dphase_gif(osol)



idx1, idx2, idx3 = find_idx(umap, "XXXBCD(t)"), find_idx(umap, "Lp(t)"),  find_idx(umap,"LpAK(t)")
plot(osol[idx1,:],osol[idx2,:],osol[idx3,:])


plot(osol[2,1:22270], osol[end,1:22270], seriesalpha = range(0,stop=1,length=22270), linewidth = range(0,10,length = 22270))

function make_testgif()
    yvals = sin.(range(0,stop=100,length=100))
    p1 = plot(1)
    anim = @animate for i in 1:100
        push!(p1, i, yvals[i])
    end
    gif(anim, "test.gif")
end

make_testgif()