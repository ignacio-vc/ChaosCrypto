#Modulo Lorenz
__precompile__(true)

module Lo

include("Automata.jl")
using AD

function generarTaylor(condIni, funcion)
    funcion(condIni)
end

function generarTaylor(condIni, funcion, t)
    funcion(condIni, t)
end

function generarSerie(polTalor)
    x = Float64[]
    push!(x, polTalor.ini)
    valores = copy(polTalor.coef)
    for i in 1:length(valores)
        push!(x, valores[i]/i)
    end
    x
end

function generaIntervalo(lista)
    p = length(lista)
    h = lista[end]
    while (lista[p] == 0)
        p -= 1
        h = lista[p]
    end
    ϵ = 1e-3
    (ϵ/h)^(1/p)
end

function horner(x,h = 1e-3)
    p = length(x)
    xt = zeros(x)
    xt[p-1] = x[p-1] + h*x[p]
    for i in 3:p
        xt[p-i+1] = x[p-i+1] + h*xt[p-i+1+1]
    end
    xt[1]
end

σ = 10
ρ = 60
β = 8/3

function lorenz(xx, t)
    x, y, z, xr, yr, zr= xx
    AD.Taylor(x, [σ])*(AD.Taylor(x,[y])-AD.Taylor(x,[x])),
    AD.Taylor(y,[x])*(AD.Taylor(y,[ρ])-AD.Taylor(y,[z]))-AD.Taylor(y,[y]),
    AD.Taylor(z,[x])*AD.Taylor(z,[y])-AD.Taylor(z,[β])*AD.Taylor(z,[z]),
    AD.Taylor(xr,[σ])*(AD.Taylor(xr,[y])-AD.Taylor(xr,[xr])),
    AD.Taylor(yr,[ρ])*AD.Taylor(yr,[xr])-AD.Taylor(yr,[yr])-AD.Taylor(yr,[xr])*AD.Taylor(yr,[zr]),
    AD.Taylor(zr,[xr])*AD.Taylor(zr,[yr])-AD.Taylor(zr,[β])*AD.Taylor(zr,[zr])
end

function integrador(x0, f, tf)
    a = generarTaylor(x0,f)
    b = map(generarSerie,a)
    suma = map(horner,b)
    sol = Array{Float64,1}[x0, [suma...]]
    t = [0.,1e-3]
    while t[end] < tf
        a = generarTaylor(sol[end],f)
        b = map(generarSerie,a)
        suma = map(horner,b)
        push!(sol,[suma...])
        push!(t,t[end]+1e-3)
    end
    sol,t
end

