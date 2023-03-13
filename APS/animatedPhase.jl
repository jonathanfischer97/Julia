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

    #previously hidden NERDSS reactions
    # (ka2,kb2), Lp + AK <--> LpAK
    # (ka2,kb2), Lp + AP <--> LpAP
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y 

function parse_input(input::AbstractString)::Vector{Float64}
    output = Float64[]
    input = replace(input, r"[\{\}]"=> "") # remove curly braces
    items = split(input, ", ")
    for s in items
        val = split(s, ": ")[2]
        push!(output, parse(Float64, val))
    end
    return output
end

input = "{'ka1': 3.6396514622137452, 'kb1': 10.311889380044056, 'kcat1': 361.5709556626289, 'ka2': 89.4797815238423, 'kb2': 500.0, 'ka3': 98.3482118725635, 'kb3': 500.0, 'ka4': 9.746019998479238, 'kb4': 121.60536126369898, 'ka7': 72.19871144097954, 'kb7': 124.66417463084926, 'kcat7': 80.47863850796283, 'y': 1500}"
p = parse_input(input)
p = [2.87928701e-01, 6.07061647e+01, 3.20457462e+02, 5.27083541e-01,
2.60234634e+00, 3.19698731e+01, 1.14855450e+02, 7.56209844e-01,
3.08201450e-01, 6.31997319e-01, 2.26216515e+01, 1.12801964e+01,
1.31157485/0.001]
pmap = (:ka1 => p[1], :kb1 => p[2], :kcat1 => p[3], :ka2 => p[4], :kb2 => p[5], :ka3 => p[6], :kb3 => p[7], :ka4 => p[8], :kb4 => p[9], :ka7 => p[10], :kb7 => p[11], :kcat7 => p[12], :y => p[13])

umap = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.9, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]#, :AK => 0., :AP => 0.]
umap = symmap_to_varmap(trimer_rn, umap)


tspan = (0.0, 100.0)

oprob = ODEProblem(trimer_rn, umap, tspan, p)
osol = solve(oprob, saveat = 0.001)#, Tsit5(), reltol=1e-8, abstol=1e-12)
plot(osol,size=(1000,600),lw=3,legend=:topleft)



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
    gif(gif1, "./lipid_timeseries.gif")
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


function find_idx(p, varmode = true)
    varmode ? map = species(trimer_rn) : map = parameters(trimer_rn) 
    for (i, k) in enumerate(map)
        if string(k) == p
            return i
        end
    end
end

for x in species(trimer_rn)
    println(x)
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
    tframevec = [Int(round(f)) for f in tframe]
    color_gradient = cgrad(:plasma, 6, categorical = true)
    pcopy = copy(p)
    idx = find_idx(param, false)
    xidx, yidx = find_idx(x), find_idx(y)

    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        pcopy[idx] = pcopy[idx]*scalar
        newprob = ODEProblem(trimer_rn, umap, tspan, pcopy)
        newsol = solve(newprob, saveat = 0.001)
        newsolarray[i] = newsol
    end

    # println(newsolarray)

    prog = Progress(length(tframevec), 1, "Making gif: ")
    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot(osol[xidx,1:f], osol[yidx,1:f], c = color_gradient[1], seriesalpha = alpharange, xlabel = "Sticky Lipid (uM)", ylabel = "Adaptor (uM)", label = "", xlims = (0,3.0), ylims = (0.5,0.9))
        scatter!((osol[xidx,f],osol[yidx,f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = "", legend = :topright)
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[xidx,1:f],nsol[yidx,1:f],c = c, seriesalpha = alpharange, label = "")
            scatter!((nsol[xidx,f],nsol[yidx,f]), markershape = :circle, markersize = 7, c=c, label = "V/A" * " * $(round(scalar^i;digits=1))")
        end
        next!(prog)
    end  
    gif(anim, "APS/"*x*"-"*y*"_"*param*"_2dphase.gif")
end

make_2dphase_gif(osol, "Lp(t)", "A(t)", "y", 0.8)




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
    gif(anim, "/APS"*x*"-"*y*"_"*param*"_"*"2dphase.gif")
end

make_2dphase_gif(osol, "XXXBCD(t)", "A(t)", "ka2", 1.2)


















function make_3dphase_gif(osol, x, y, z, param, scalar) #make gif of #D phase plot 
    tframe = range(2,stop=length(osol.t),length=500)
    tframevec = [Int(round(f)) for f in tframe]
    color_gradient = cgrad(:plasma, 6, categorical = true)
    pcopy = copy(p)
    idx = find_idx(param, false)
    xidx, yidx, zidx = find_idx(x), find_idx(y), find_idx(z)


    newsolarray = Array{Any}(undef,length(color_gradient)-1)
    @views for i = 1:length(color_gradient)-1
        pcopy[idx] = pcopy[idx]*scalar
        newprob = ODEProblem(trimer_rn, umap, tspan, pcopy)
        newsol = solve(newprob, saveat = 0.001)
        newsolarray[i] = newsol
    end

    prog = Progress(length(tframevec), 1, "Making gif: ")
    anim = @animate for f in tframevec
        alpharange = collect(((10^x)/10 for x in range(0,stop=1,length=f)))
        # widthrange = collect(((10^x)/10 for x in range(0,stop=5,length=f)))

        # plot(osol[2,1:f],osol[end,1:f], seriesalpha = alpharange, linewidth = widthrange, label = "V/A", c = color_gradient[1])#,xlims=(0,maximum(osol[2,:])*1.1),ylims=(0,maximum(osol[end,:])*1.1), legend = :topleft)
        plot3d(osol[xidx,1:f], osol[yidx,1:f], osol[zidx,1:f], arrow = true, c = color_gradient[1], seriesalpha = alpharange)#, xlims=(0.12,0.2),ylims=(0,0.22), label = param, camera=(360k/100, 30))
        # scatter3d!((osol[2,f],osol[end,f],osol[4,1:f]), markershape = :circle, markersize = 7, color = color_gradient[1], label = "V/A")
        @views for (i,c) in enumerate(color_gradient[2:end])
            nsol = newsolarray[i]

            plot!(nsol[xidx,1:f],nsol[yidx,1:f],nsol[zidx,1:f], arrow = true, c = c, seriesalpha = alpharange, label = param * " * $(round(scalar^i;digits=1))")
            # scatter3d!((nsol[2,f],nsol[end,f],nsol[4,1:f]), markershape = :circle, markersize = 7, c=c, label = "V/A * $(round(1.1^i;digits=1))")

        end
        next!(prog)
    end  
    gif(anim, "APS/"*x*"_"*y*"_"*z*"_3dphase.gif")
end

make_3dphase_gif(osol, "Lp(t)", "A(t)", "LpAK(t)", "y", 0.8)


