using Plots
using LinearAlgebra

# 1. Definir la curva C (circunferencia en el plano s)
N = 1000
radius = 10
center = 12
angles = LinRange(0, 2π, N+1)[1:end-1]  # para evitar repetir el punto inicial
s = center .+ radius * exp.(-1im .* angles)

# 2. Definir la función L(s)
L(s) = (2s^2 + 4s + 6.5) / (s^2 + 5.5s + 2.5)

# 3. Calcular la imagen de la curva
Lc = L.(s)
L1c = 1 .+ Lc

# 4. Función para contar las vueltas alrededor de un punto
function winding_number(curva, punto)
    arg_diffs = diff(angle.(curva .- punto))
    # Corregimos saltos de ±2π
    arg_diffs = map(x -> mod(x + π, 2π) - π, arg_diffs)
    return round(sum(arg_diffs) / (2π))
end

# 5. Calcular número de vueltas
vueltas_0 = winding_number(L1c, 0)     # 1 + L(C) alrededor de 0
vueltas_m1 = winding_number(Lc, -1)    # L(C) alrededor de -1

# 6. Graficar ambas curvas
plot(Lc, label="L(C)", color=:blue, aspect_ratio=:equal)
plot!(L1c, label="1 + L(C)", color=:red)
scatter!([0], [0], label="origen", color=:black, marker=:x)
scatter!([-1], [0], label="-1", color=:green, marker=:x)
title!("Curvas L(C) y 1 + L(C)")
xlabel!("Re")
ylabel!("Im")

println("Vueltas de 1 + L(C) alrededor de 0: $vueltas_0")
println("Vueltas de L(C) alrededor de -1: $vueltas_m1")
