using Statistics, Dates, Random
using Plots, LaTeXStrings, BenchmarkTools

function Heun(x, dt, dx, T)
    rng = MersenneTwister(12)
    noise = randn!(rng, zeros(length(x))) * sqrt(2 * T / (dt * dx))
    yt = x + f(x, dx) * dt + noise * dt
    y = x + (f(x, dx) + f(yt, dx)) * dt / 2 + noise * dt
    return y
end

function f(x, dx)
    y = (circshift(x, 1) + circshift(x, -1) - 2 * x) / (dx)^2
    return y
end

function sqv(vec)
    return vec .^ 2
end

function time_evoluition(h, t, M_measure, dt, dx, T)
    order_final = zeros(length(t))
    order_mean = zeros(length(t))
    order2mean = zeros(length(t))
    for n ∈ eachindex(M_measure)
        for i ∈ 1:length(t)-1
            h_pre = h[i]
            h_next = h[i+1]
            h_next .= Heun(h_pre, dt, dx, T)
            # order_mean[i] = mean(h_pre)
            # order2mean[i] = mean(sqv.(h_pre))
        end
        # order_mean[end] = mean(h[end])
        # order2mean[end] = mean(sqv.(h[end]))
        order_mean = [mean(h[j]) for j ∈ 1:length(t)]
        order2mean = [mean(sqv.(h[j])) for j ∈ 1:length(t)]
        order = sqrt.(order2mean - sqv.(order_mean))
        order_final = order_final + order
    end
    println(h[end])
    order_final = order_final / M_measure
    return order_final
end

# parameters
L = 10
dx = 1e-1
x = 0:dx:L

t_max = 100
dt = 1e-3
t = 0:dt:t_max

T = 1

M_measure = 1e3

# initial condition
h0 = zeros(length(x))

h = [zeros(length(x)) for i ∈ t]

# main
order_final = time_evoluition(h, t, M_measure, dt, dx, T)

plot(t, order_final,
    framestyle=:box,
    xlabel=L"t",
    ylabel=L"order-parameter",
    # xlims = [0, 75],
    label=L"T=1",
    #ls = :dashdot,
    lw=2,
    c=2,
    #la = 0.5,
    #xtick = ([0, 50], [1, 5]),
    title=L"growth"
)

plot(log.(t), log.(order_final),
    framestyle=:box,
    xlabel=L"t",
    ylabel=L"order-parameter",
    # xlims = [0, 75],
    label=L"T=1",
    #ls = :dashdot,
    lw=2,
    c=2,
    #la = 0.5,
    #xtick = ([0, 50], [1, 5]),
    title=L"growth"
)