using Plots 
using Catalyst
using DifferentialEquations
using Symbolics
using Latexify

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





input = "{'ka1': 3.6396514622137452, 'kb1': 10.311889380044056, 'kcat1': 361.5709556626289, 'ka2': 89.4797815238423, 'kb2': 500.0, 'ka3': 98.3482118725635, 'kb3': 500.0, 'ka4': 9.746019998479238, 'kb4': 121.60536126369898, 'ka7': 72.19871144097954, 'kb7': 124.66417463084926, 'kcat7': 80.47863850796283, 'VA': 1500}"
p = [1.00000000e+01,  1.36182339e+02,  3.20457462e+02,  5.27083541e-01,
2.60234634e+00,  3.19698731e+01,  1.14855450e+02,  7.56209844e-01,
3.08201450e-01,  6.31997319e-01,  2.26216515e+01,  1.36182339e+01,
1.31157485/0.001]

p = parse_input(input)

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.9, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]

tspan = (0.0, 100.0)

oprob = ODEProblem(trimer_rn, u0, tspan, p)
osol = solve(oprob, Tsit5(), reltol=1e-8, abstol=1e-12)
plot(osol,size=(1000,600),lw=3,legend=:topleft)