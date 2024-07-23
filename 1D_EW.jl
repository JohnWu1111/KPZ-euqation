using Statistics, Dates
using Plots, LaTeXStrings,BenchmarkTools

# parameters
L = 10
dx = 5e-2
x = 0:dx:L

t_max = 10
dt = 1e-2
t = 1:dt:t_max

# initial condition
h0 = zeros(length(x))

h = [zeros(length(x)) for i ∈ t]
