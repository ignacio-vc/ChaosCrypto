using ChaosCrypto
using PyPlot
function lorenz(xx)
    x, y, z, xr, yr, zr= xx
    Taylor(x, [σ])*(Taylor(x,[y])-Taylor(x,[x])),
    Taylor(y,[x])*(Taylor(y,[ρ])-Taylor(y,[z]))-Taylor(y,[y]),
    Taylor(z,[x])*Taylor(z,[y])-Taylor(z,[β])*Taylor(z,[z]),
    Taylor(xr,[σ])*(Taylor(xr,[y])-Taylor(xr,[xr])),
    Taylor(yr,[ρ])*Taylor(yr,[xr])-Taylor(yr,[yr])-Taylor(yr,[xr])*Taylor(yr,[zr]),
    Taylor(zr,[xr])*Taylor(zr,[yr])-Taylor(zr,[β])*Taylor(zr,[zr])
end
vecs,t = integrador([1.,1.,1.,1.,100.,100.],lorenz,10.)
title("Ejercicio 1.1")
plot(t, [x[2] for x in vecs], ".", label=L"$y$")
plot(t, [x[5] for x in vecs], label=L"$y_r$")
legend()
