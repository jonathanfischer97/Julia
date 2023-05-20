"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(p,  prob::ODEProblem)
    Y = nothing
    try 
        Y = solve(remake(prob, p=p), saveat=0.1, save_idxs=1, maxiters=10000, verbose=false)
        if Y.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) || any(x==1 for array in isnan.(Y) for x in array) || any(x==1 for array in isless.(Y, 0.0) for x in array)
            return 1.0
        end
    catch e 
        if e isa DomainError #catch domain errors
            return 1.0
        else
            rethrow(e) #rethrow other errors
        end
    end
    fitness, period, amplitude = CostFunction(Y)


    return (fit = -fitness, per = period, amp = amplitude)
end