include("RK4.jl")
using PyPlot
using RK
σ=16.
ρ=45.6
β=4.

function lorenzSin(xx,t)
    x, y, z, xr, yr, zr = xx
    
    φ = 1e-3
    m = sin(φ*t)
    s = x + m
    
    [σ*(y-x), (ρ*x - y - x*z), x*y-β*z,
     σ*(yr-xr), (ρ*s - yr - s*zr) , (s*yr - β*zr)]
end

xs, ts = integrar(lorenzSin,[1., 1., 1., 10., 10., 10.],0. ,100. ,1e-3);
m=sin(ts)
s=[x[1] for x in xs]+m
mhat = s-[x[4] for x in xs]
plot(ts,m)
plot(ts,mhat)
