__precompile__(true)

module ChaosCrypto

export AD, LO

export Taylor, paso2,paso1
export igualdad,logo,expo,seno,coseno
export generarTaylor,generarTaylor,generarSerie
export generaIntervalo,horner,integrador
export σ,ρ,β

include("Automata.jl")
include("Lorenz.jl")

end
