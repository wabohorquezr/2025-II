using ControlSystems
using RobustAndOptimalControl
using LinearAlgebra
using Plots
#
#Definción de la función L
L = zpk([],[-1, -2],2)
fig = nyquistplot(L; unit_circle=true)  
re_n, im_n, w = nyquistv(L)  
fig_n = nyquistplot(L,xlims=(-1.5, 1), ylims=(-1.5, 1.5),linecolor="blue")
fig_2 = nyquistplot!(L,-w,xlims=(-1, 1.5), ylims=(-1.5, 1.5),linecolor="blue")
#
#Definción de la función L
L = zpk([],[-1, -2],2)
fig = nyquistplot(L; unit_circle=true)  
re_n, im_n, w = nyquistv(L)  
fig_n = nyquistplot(L,xlims=(-1.5, 1), ylims=(-1.5, 1.5),linecolor="blue")
fig_2 = nyquistplot!(L,-w,xlims=(-1, 1.5), ylims=(-1.5, 1.5),linecolor="blue")
#




L1 = zpk([],[-0.5,-1, -2],1)
fig1 = nyquistplot(L1; unit_circle=true)  
re_n, im_n, w = nyquistv(L1) 
w_n = [reverse(-w);w] 
fig_n = nyquistplot(L1,xlims=(-1.5, 1), ylims=(-2, 2),linecolor="blue")
fig_2 = nyquistplot!(L1,-w,xlims=(-1.5, 1), ylims=(-2, 2),linecolor="blue")
#
fig_3 = nyquistplot!(11.25*L1,xlims=(-3, 1), ylims=(-2, 2),linecolor="red")
fig_4 = nyquistplot!(11.25*L1,-w,xlims=(-3, 1), ylims=(-2, 2),linecolor="red")
#
# Estabilidad interna
P = zpk([],[0,2],1)
K = zpk([2],[-2],1)
S, PS, CS, T = gangoffour(P, K; minimal=true)
nyquistplot(P*K,xlims=(-1.5, 1), ylims=(-2, 2))
L=zpk([],[0, -1, -2],1)
isstable(S)
isstable(PS)
isstable(CS)
isstable(T)

# Retardos
L = zpk([],[-1,-2,-2.5],5)

margin(L)
nyquistplot(L)
nyquistplot!(9.45L)
nyquistplot!(5L)
k = 1
T = feedback(k*L,1)
plot(step(T))
for k in [5, 9.4]
end
Lk = TransferFunction[k*L for k in [1, 5, 9.4]]
T = feedback.(Lk,1)
plot(step.(T))

nyquistplot(L)
nyquistplot!(L*delay(0.5))
nyquistplot!(L*delay(1))
nyquistplot!(2L*delay(1))

plot(step(feedback(L,1)))
plot!(step(feedback(L*delay(0.5),1)))
plot!(step(feedback(L*delay(1),1)))
plot!(step(feedback(2*L*delay(1),1)))

T_delay = T * delay(1.2)
marginplot(5T_delay)
plot(step(T_delay))
plot!(step(T))
nyquistplot(T_delay)
S, PS, CS, T, RY, RU, RE = gangofseven(P,K,tf([1],[1]))
sys_P =[P P 1]
sys_K = [1 1;P P]
sys_lc = feedback(sys_P,sys_K;U1=:[3],W1=:[1 2],U2=:[1],W2=:[2])
fig = nyquistplot(L, xlims=(-1.5, 1), ylims=(-2, 1); unit_circle=true)  
plot(fig)

#
marginplot(L)


A = [1 2; 3 4];
C = [5 6; 7 8];
X = lyap(A, C)
eigen(A)
eigen(X)
using MatrixEquations
X = lyapc(A,C)
eigen(X)
#
A = [-1 2;-2 -1]
eigen(A)
Q = I(2)
X = lyapc(A,Q)
eigen(X)
X = lyap(A,Q)

using ControlSystemsBase, Plots
G = tf(1, [1, 1, 1])
res = step(G, 20) # Simulate 20 seconds step response
plot(res)

si = stepinfo(res)
plot(si)
diagm(0 => [1, 0, 2])
rank(diagm(0 => [1, 0, 2]))

using ControlSystems
using LinearAlgebra: I
using Plots

A = [0 1; 0 0]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0)
Q = I
R = I
L = lqr(sys,Q,R)

u(x,t) = -L*x .+ (t > 4) # State feedback + step disturbance
t  = 0:0.1:12
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x', lab=["Position" "Velocity"], xlabel="Time [s]"); vline!([4], lab="Step disturbance", l=(:black, :dash, 0.5))

#
using ControlSystems
s = tf('s')
P = zpk([],[-1, -0.25], 1) 
k = 1
K = k*(1 + 1/(8s))
L = K * P 
S = 1/(1 + L)
minreal(S)
T = L/(1+L)
pole(L)
pole(T)
tzero(L)
tzero(T)

#
k = 5
K = k*(1 + 1/(8s))
L = K * P 
S = 1/(1 + L)
minreal(S)
T = feedback(L,1)
pole(L)
pole(T)
tzero(L)
tzero(T)

#
plot(step(T))

plot(step(S))

include("ControlUN.jl")
step_plt(T,12.5)