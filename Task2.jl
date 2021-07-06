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
function CA(init_state, max_time, stop_time, p_gamma, p_xi)
    t = max_time + 1
    A = zeros(t, 2 * t + 1)

    if init_state == 0
        A[1, t+1] = 1
    else
        A[1, 2:2*max_time+2] = init_state
    end
    
    for time = 1:stop_time
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
#the output showed activity in the range [0.6,0.8]
#expected_endstate was then run for each p in the range [0.6:0.01:0.8]
function expected_endstate(samp, p)
    v = 0
    for i = 1:samp
        v += CA(0, 2000, 2000, 0.9, p)[2001, 2002]
    end
    return v / samp
end

#normalized hamming function
function ham(max_time, stop_time, p, init_dif)
    e = zeros(2 * max_time + 1)
    while e'e < init_dif
        e[rand(1:1:length(e))] = 1
    end

    v = rand(Bernoulli(0.5), 2 * max_time + 1)
    vc = copy(v)
    for i = 1:length(e)
        vv[i] = (vv[i] - e[i])^2
    end

    A0 = CA(v, max_time, stop_time, p, 0.9)
    A1 = CA(vv, max_time, stop_time, p, 0.9)
    m, n = size(A0)

    D = A0 - A1
    H = zeros(stop_time)
    for i = 1:stop_time
        H[i] = n - count(x -> (x == 0), D[i, :])
    end

    return H / (n - 2)
end
                    
#Average hamming function                   
function smooth_ham(max_time, stop_time, p, init_dif, samp)
    H = ham(max_time, stop_time, p, init_dif)
    for j = 1:samp-1
        H += ham(max_time, stop_time, p, init_dif)
    end
    return H / samp
end

