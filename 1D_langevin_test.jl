using Statistics, Dates
using Plots, LaTeXStrings,BenchmarkTools

function Heun(x, dt, T)
    noise = randn() * sqrt(2*T / dt)
    y1t = x[1] + x[2] * dt
    y2t = x[2] + (-x[1] - x[2] + noise) * dt
    y = zeros(2)
    y[1] = x[1] + (x[2] + y2t) * dt / 2
    y[2] = x[2] + (-x[1] - x[2] - y1t - y2t + 2 * noise) * dt / 2
    return y
end

function sqv(vec)
    return vec.^2
end

#@time begin

starttime = now()

q0 = [0; 1]
dt = 1e-2
t_max = 100
t = 0:dt:t_max
nt = length(t)
t_max = 100
T = 1:t_max
nT = length(T)
N_total = 100

q_mean = [zeros(2) for i = 1:nT]

qtv = [zeros(2) for i = 1:nt]
first(qtv) = q0

for m ∈ 1:N_total

for n ∈ eachindex(T)   
    for i in 1:nt-1
        qt = qtv[i]
        qtp = qtv[i+1]
        qtp .= Heun(qt, dt, T[n])
    end
    q_mean[n] = q_mean[n] + mean(sqv.(qtv[round(Int,nt/10):end]))
    # for i in 1:nt-1
    #     qt[:, i+1] .= Heun(qt[:, i], dt, T[n])
    # end
end
end
q_mean = q_mean./N_total

endtime = now()
duration = endtime - starttime

#aa = 1
plot(T,[q_mean[i][1] for i ∈ eachindex(q_mean)],
framestyle = :box,
xlabel = L"T",
ylabel = L"E",
# xlims = [0, 75],
label = L"V",
#ls = :dashdot,
lw = 2,
c = 2,
#la = 0.5,
#xtick = ([0, 50], [1, 5]),
title = L"growth"
)
plot!(T,[q_mean[i][2] for i ∈ eachindex(q_mean)],
label = L"T",
)

#end
