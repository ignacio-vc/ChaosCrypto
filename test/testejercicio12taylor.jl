using ChaosCrypto
using PyPlot
function lorenz12(xx,t)
    x, y, z, xr, yr, zr= xx
    AD.Taylor(x, [LO.σ])*(AD.Taylor(x,[y])-AD.Taylor(x,[x])),
    AD.Taylor(y,[x])*(AD.Taylor(y,[LO.ρ])-AD.Taylor(y,[z]))-AD.Taylor(y,[y]),
    AD.Taylor(z,[x])*AD.Taylor(z,[y])-AD.Taylor(z,[LO.β])*AD.Taylor(z,[z]),
    AD.Taylor(xr,[LO.σ])*(AD.Taylor(xr,[y])-AD.Taylor(xr,[xr])),
    AD.Taylor(yr,[LO.ρ])*AD.Taylor(yr,[xr])-AD.Taylor(yr,[yr])-AD.Taylor(yr,[xr])*AD.Taylor(yr,[zr]),
    AD.Taylor(zr,[xr])*AD.Taylor(zr,[yr])-AD.Taylor(zr,[LO.β])*AD.Taylor(zr,[zr])
end
vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenz12,10.)
title("Ejercicio 1.2")
plot3D([x[1] for x in vecs],[x[2] for x in vecs],[x[3] for x in vecs])
legend()
