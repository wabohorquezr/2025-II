include("ControlUN.jl")
using LinearAlgebra
A = [-6 -3.5;
      6  4]
B=[0; 1]
# Localización deseada
polos=[-3, -1]
# Calcula matriz de ganancia para localizar los polos
G=place(A,B,polos)
ALC = A - B*G
eigen(ALC)
