using Plots 
using Catalyst
using DifferentialEquations
using Symbolics
using Latexify
using Statistics


function calc_D(product, Ddict)
    # Filter out the letter "p" from the product string and split the remaining characters
    substrings = filter(c -> c ≠ 'p', split(product, ""))

    # Initialize a list to store the average diffusion constants
    avg_Ds = []

    # Loop over the substrings and compute the average diffusion constant for each one
    for substring in substrings
        # Use the substring as a key to access the diffusion constant value in the Ddict dictionary
        Dvals = Ddict[substring]

        # Compute the average diffusion constant over the x, y, and z directions
        avg_D = (Dvals[1] + Dvals[2] + Dvals[3]) / 3

        # Append the average diffusion constant to the list
        push!(avg_Ds, avg_D)
    end

    # Compute the average relative diffusion constant over all the substrings
    avg_rel_D = mean(avg_Ds)

    # Return the result
    return avg_rel_D
end




# Define the diffusion constant dictionary
Ddict = Dict("L" => [0.5, 0.5, 0.0],
             "K" => [25.0, 25.0, 25.0],
             "P" => [25.0, 25.0, 25.0],
             "A" => [25.0, 25.0, 25.0])

# Calculate the average relative diffusion constant for the "AKL" product
avg_rel_D = calc_D("L", Ddict)

# Print the result
println("Average relative diffusion constant for AKL: ", avg_rel_D)




trimer_rn = @reaction_network trirn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1y*y,kb1y), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7y*y,kb7y), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 

    #previously hidden NERDSS reactions
    (ka2i,kb2i), Lp + AK <--> LpAK
    (ka2y*y,kb2y), Lp + AKL <--> LpAKL
    (ka2i,kb2i), Lp + AP <--> LpAP
    (ka2y*y,kb2y), Lp + APLp <--> LpAPLp
    (ka3i,kb3i), A + K <--> AK
    (ka4i,kb4i), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3y*y,kb3y), LpA + LK <--> LpAKL
    (ka4y*y,kb4y), LpA + LpP <--> LpAPLp
    (ka1i,kb1i), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7i,kb7i), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y ka1i kb1i ka2i kb2i ka3i kb3i ka4i kb4i ka7i kb7i ka1y kb1y ka2y kb2y ka3y kb3y ka4y kb4y ka7y kb7y 

osys = convert(ODESystem, trimer_rn)

for eq in osys.eqs
    println(eq)
end

function format_equations(eqs::Symbolics.Arr{Equation,1})
    for ode in eqs
        ode_str = string(ode)
        # ode_str = replace(ode_str, "~" => "=")
        ode_str = replace(ode_str, "(t)" => "")
        ode_str = replace(ode_str, "Differential(" => "d")
        ode_str = replace(ode_str, ") ~" => " =")
        println(ode_str)
    end
end





eqs = format_equations(osys.eqs)


p = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.]

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.6, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]