using ChaosCrypto
using PyPlot

function lorenzSin(xx,t)
    x, y, z, xr, yr, zr = xx
    
    φ = 1e-3
    m = sin(φ*t)
    s = x + m
    
    [Taylor(x,[σ*(y-x)]), 
        Taylor(y,[(ρ*x - y - x*z)]), 
        Taylor(z,[x*y-β*z]),
        Taylor(xr,[σ*(yr-xr)]), 
        Taylor(yr,[(ρ*s - yr - s*zr)]), 
        Taylor(zr,[(s*yr - β*zr)])]
end

vecs,t = integrador([1.,1.,1.,1.,100.,100.],lorenzSin,10.)

m = sin(t)
s = [x[1] for x in vecs] + m
mhat = s-[x[4] for x in vecs]
plot(t,m)
plot(t,mhat)
