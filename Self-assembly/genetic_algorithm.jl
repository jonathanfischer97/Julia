













## COST FUNCTIONS and dependencies 
function getDif(indexes, arrayData)
    arrLen = length(indexes)
    println(1)
    sum = 0
    for (i, ind) in enumerate(indexes)
        println(i)
        if i == arrLen
            println(3)
            break 
        end
        sum += arrayData[ind] - arrayData[indexes[i+1]]
    end
    sum += arrayData[indexes[end]] 
    return sum #return statement is unneeded but just covering bases 
end
    
function getSTD(indexes, arrayData, window) #get standard deviation of fft peak indexes
    numPeaks = length(indexes)
    arrLen = length(arrayData)
    sum = 0
    for ind in indexes 
        minInd = max(1, ind - window)
        maxInd = min(arrLen, ind+window)
        sum += std(arrayData[minInd:maxInd])
    end
    sum = sum/numPeaks 
    return sum
end
 
function getFrequencies1(y) #normally would slice y by interval for a sampling rate but not here
    y1 = y
    #fft sample rate: 1 sample per 5 minutes
    res = broadcast(abs,rfft(y1)) #broadcast abs to all elements of rfft return array. Think .= syntax does the same 
    #normalize the amplitudes
    norm_res = res/cld(1000, 2)
    return norm_res #smallest integer larger than or equal to. Rounding up
end

##This cost function is the goal but doesn't yet work, IGNORE
function costTwo(y)
    Y = Array(solve(prob, tspan = (0., 100.), p = y, Tsit5()))
    p1 = Y[1,:]
    fftData = getFrequencies1(p1)

    indexes, _ = findmaxima(fftData) 
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return std + diff
end

#homemade peakfinder
function findlocalmaxima(signal::Vector)
    inds = Int[]
    buff = Zygote.Buffer(inds) #Zygote.Buffer creates mutable array that is invisible to autodifferention, aka no matrix operations allowed
    if length(signal)>1
        if signal[1]>signal[2]
            push!(buff,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(buff,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(buff,length(signal))
        end
    end
    return copy(buff) #return copy of mutable buffer array to make immutable (essential for Zygote pullback)
  end

#current testing cost function, USE 
function testcost(p)
    Y = Array(solve(prob, tspan = (0., 100.), p = p, Tsit5()))
    p1 = Y[1,:]
    fftData = getFrequencies1(p1)
    indexes = findlocalmaxima(fftData)[1]
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return std + diff
end
