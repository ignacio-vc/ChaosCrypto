include("testejercicio2.jl")
title("Ejercicio 9.6.4 - 2, X vs Xr")
plot(t22, [x[1] for x in vecs22], label="x")
plot(t22, [x[4] for x in vecs22], label="xr")
legend()
show()