#Ejercicio 9.6.5

using ChaosCrypto
using PyPlot

<<<<<<< HEAD
fig = figure()

@manipulate for ρr in [1,10,24,25,60]
    withfig(fig) do
        function lorenzSin(xx,t)
            x, y, z, xr, yr, zr = xx
            
            φ = 1e-3
            m = sin(φ*t)
            s = x + m
            
            [AD.Taylor(x,[LO.σ*(y-x)]), 
                AD.Taylor(y,[(ρr*x - y - x*z)]), 
                AD.Taylor(z,[x*y-LO.β*z]),
                AD.Taylor(xr,[LO.σ*(yr-xr)]), 
                AD.Taylor(yr,[(ρr*s - yr - s*zr)]), 
                AD.Taylor(zr,[(s*yr - LO.β*zr)])]
        end
        vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenzSin,10.)

        m = sin(t)
        s = [x[1] for x in vecs] + m
        mhat = s-[x[4] for x in vecs]
        title("Comparación señal y resultado")
        plot(t,m, label="Señal")
        plot(t,mhat, label="Resultado")
        legend()
    end
end
=======
function lorenzSin(xx,t)
    x, y, z, xr, yr, zr = xx
    
    # x,y,z gobernados por: ρ = 60, σ = 10, β = 8/3 (default)
    # xr,yr,zr gobernados por parametros de ejercicio 9.6.2
    # Comenzamos a integrar numericamente

    φ = 1e-3
    m = sin(φ*t)
    s = x + m
    
    [AD.Taylor(x,[LO.σ*(y-x)]), 
        AD.Taylor(y,[(LO.ρ*x - y - x*z)]), 
        AD.Taylor(z,[x*y-LO.β*z]),
        AD.Taylor(xr,[LO.σ*(yr-xr)]), 
        AD.Taylor(yr,[(LO.ρ*s - yr - s*zr)]), 
        AD.Taylor(zr,[(s*yr - LO.β*zr)])]
end

# Escogiendo condiciones iniciales diferentes para y, yr; z, zr..
vecs,t = LO.integrador([1.,1.,1.,1.,100.,100.],lorenzSin,10.)

m = sin(t)
s = [x[1] for x in vecs] + m
mhat = s-[x[4] for x in vecs]
title("Ejercicio 3 con Taylor")
plot(t,m)
plot(t,mhat)
show()
>>>>>>> sonidoCC
