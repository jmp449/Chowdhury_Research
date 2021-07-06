using Plots
using StatsBase
using JLD
using Distributions
using DelimitedFiles
using LaTeXStrings
using Colors
using CurveFit
using GR

#--------------------------
# heaviside step function
function H(x)
    if x >= 1
        return 1
    else
        return 0
    end
end

#produces state space of cellular autamaton for a given rule
function CA(inn, max_time, stop_time, p_gamma, p_xi)
    t = max_time + 1
    A = zeros(Int64, (t, 2 * t + 1))

    if inn == 0
        A[1, t+1] = 1
    else
        A[1, 2:2*max_time+2] = inn
    end

    for time = 1:stop_time
        for pos = -time:time
            #C(r,t)
            c00 = A[time, pos+1+t]
            #C(r+1,t)
            cp1 = A[time, pos+2+t]
            #C(r-1,t)
            cn1 = A[time, pos+t]
            #C(r,t+1)
            A[time+1, pos+1+t] =
                c00 * (1 - rand(Bernoulli(p_xi))) +
                (rand(Bernoulli(p_gamma))) * (1 - c00) * H(cp1 + cn1)
        end
    end
    return A
end
#returns sample of final states of r=0 cell at t=2000 for samp=200 and a given p
#expected_endstate was run for each p in the range [0:0.01:1]
function expected_endstate(samp, p)
    v = 0
    for i = 1:samp
        v += CA(0, 2000, 2000, p, 0.9)[2001, 2002]
    end
    return v / samp
end
#shift deals with nonegative indexing in julia
function shift(M, x)
    return M[convert(Int64, floor(10 * x + 1))]
end

function average(x, num)
    D = copy(x)
    for i = 1:num
        for j = 2:length(x)-1
            D[j] = 0.5 * (D[j+1] + D[j])
        end
    end
    return D
end

#normalized hamming function
function ham(max_time, stop_time, p, init_dif)
    e = zeros(2 * max_time + 1)
    while e'e < init_dif
        e[rand(1:1:length(e))] = 1
    end

    v = rand(Bernoulli(0.5), 2 * max_time + 1)
    vv = copy(v)
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

    return H / n
end
function smooth_ham(max_time, stop_time, p, init_dif, samp)
    H = ham(max_time, stop_time, p, init_dif)
    for j = 1:samp-1
        H += ham(max_time, stop_time, p, init_dif)
    end
    return H / samp
end
#-------------------------------------
#Average value plots below
v = zeros(11)
vv = zeros(11)
for i = 0:10
    vv[i+1] = expected_endstate(100, (690 + i) / 1000)
end
plot(vv)
C(y) = shift(v, y)
plot(v)
plot(
    #xaxis = (:log10, [10^-0.02, 1]),
    #yaxis = (:log10, [10^-0.5, 1]),
    xticks = 0:0.1:1,
    xlabel = "p",
    ylabel = L"<C(r=0,t=2000)> (p)",
    C,
    legend = :none,
)
savefig("avg_linear_plot.png")
#----------------------------------
#Hamming plots below
#critical probability
pcc = 0.696
max_time = 2000
stop_time = 2000

H_l = smooth_ham(max_time, stop_time, 0.68, 100, 100)
H_eq = smooth_ham(max_time, stop_time, pcc, 100, 100)
H_g = smooth_ham(max_time, stop_time, 0.71, 100, 100)
plot(
    xlims = (1, 2000),
    ylims = (0, 0.0175),
    xaxis = (:log10, [10^0.5, stop_time]),
    yaxis = (:log10, [10^-4, 10^-1.5]),
    xlabel = L"t",
    ylabel = L"H(t)",
    legend = :topleft,
    label = L"p<p_{crit}",
    H_l,
)
plot!(H_eq, label = L"p=p_{crit}")
plot!(H_g, label = L"p>p_{crit}")
annotate!(700, 0.010, text(L"p = 0.71", :black, :right, 8))
annotate!(1900, 0.0043, text(L"p_{crit} = 0.696", :black, :right, 8))
annotate!(800, 10^-4, text(L"p = 0.68", :black, :right, 8))
savefig("ham_loglog.png")

A = CA(0, 50, 50, 1, 1)
