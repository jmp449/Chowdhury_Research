using Plots
using StatsBase
using Distributions
using DelimitedFiles

# heaviside step function
function H(x)
    if x >= 1
        return 1
    else
        return 0
    end
end

#produces matrix containing evolution of cellular autamaton
#Note: C(r,t) = A[t + 1,r + (max_time + 2)], init_state = C(r,t=0)
function CA(max_time, p_gamma, p_xi)
    t = max_time + 1
    
    A = zeros(t, 2 * t + 1)
    A[1, t+1] = 1
    
    for time = 1:max_time
        for pos = -time:time
            #C(r,t)
            c00 = A[time, pos+1+t]
            
            #C(r + 1,t)
            cp1 = A[time, pos+2+t]
            
            #C(r - 1,t)
            cn1 = A[time, pos+t]
            
            #C(r,t + 1)
            A[time+1, pos+1+t] =
                c00 * (1 - rand(Bernoulli(p_gamma))) +
                (rand(Bernoulli(p_xi))) * (1 - c00) * H(cp1 + cn1)
        end
    end
    return A
end

#returns sample of final states of r=0 cell at t=2000 for samp=200 and a given p_xi
#expected_endstate was run for each p in the range [0:0.01:1]
#the output showed activity in the range [0.6,1.0]
#expected_endstate was then run for each p in the range [0.6:0.01:1.0]
function expected_endstate(samp, p)
    v = 0
    for i = 1:samp
        v += CA(2000, 0.9, p)[2001, 2002]
    end
    return v / samp
end
