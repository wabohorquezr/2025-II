using Alert
using Plots; gr(linewidth = 2, grid=:true)
using Printf
Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)
using LaTeXStrings
#using Polynomials
using InvertedIndices
using Latexify
using PolynomialRoots
using ControlSystems
using LinearAlgebra
function dise_2p(P,T,polos_obs)
   m = length(polos_obs)
   P = tf(P);
   N, D = real.(tfdata(P));
   n = length(D) - 1
   if m < n-1
      throw("dise_2p: número de polos del observador, insuficientes")
      return
   end

     # CALCULA UN CONTROLADOR DE DOS PARÁMETROS
     # Para la planta dada por una función de transferencia P
     # La función de transferencia de lazo cerrado deseada es T
     # m: orden del controlador
     # PolosObs: Polos del observador
    
     # Forma polinomio del observador
   polos = copy(polos_obs);
   Dhb = real.(poly(polos));
   divisor = tf(1,N);
   H = minreal(T * divisor); 
   Nh, Dh = real.(tfdata(H));
   L = convol(Nh, Dhb);
   F = convol(Dh, Dhb);
   n = length(D)-1;
   m = length(polos); # El número de polos del observador
   # Forma polinomio con los polos
   l = n - length(N) + 1;
   Na = zeros(length(D));
   Na[l+1:end] = N;
   Sm = zeros(Float64,n+m+1,2*(m+1));
   Smj = zeros(Float64,n+m+1,2*m+1);
   
   for k=1:m+1
      Sm[k:k+n,k]=D;
      Sm[k:k+n,k+m+1]=Na;
   end
   #
   if m > n - 1
      # Si es posible, disena un controlador que rechaza perturbaciones 
      Smj=Sm[:,Not(m+1)]
      if rank(Smj)>=n+m+1
         XY=Smj\F;
         A_c=[XY[1:m];0.];
         M_c=XY[m+1:end];
      else
         XY=Sm\F;
         A_c=XY[1:m+1];
         M_c=XY[m+2:end];  
      end
   elseif m==n-1
      XY=Sm\F;
      A_c=XY[1:m+1];
      M_c=XY[m+2:2*m+2]; 
   end
   Kr = tf(real.(L[:,1]),real.(A_c[:,1]))
   Ky = tf(real.(M_c),real.(A_c[:,1]));
 #  T = minreal(Kr .* feedback(P, Ky));
 #  Gur=minreal(T/P);
   Kr, Ky
end

function asigne_polos(P,polos)
    # CALCULA UN CONTROLADOR CON REALIMENTACIÓN UNITARIA
    # Para la planta dada por una funcion de transferencia P
    # Los polos se asignan en el vector pol
    # m: orden del controlador
    P1 = tf(P);
    # Obtiene los datos de la planta
    P1 = tf(P);
    # Obtiene los datos de la planta
    N,D = tfdata(P1);
    # n =  orden de la planta
    n = length(D)-1;
    l = length(D)-length(N);
    # Calcula orden del sistema de lazo cerrado
    nMm = length(polos);
    m = nMm-n;
    # Forma polinomio con los polos
    DT = poly(polos);
    Na = zeros(length(D));
    Na[l+1:end]  =  N;
    # Forma la Matriz de Sylvester
    Sm = zeros(n+m+1,2*(m+1));
    for k = 1:m+1
      Sm[k:k+n,k] = D;
      Sm[k:k+n,k+m+1] = Na;
   end
   ind_error = 0;
   if m > n - 1;
      # Si es posible, disena un controlador que rechaza perturbaciones 
      #Smj = Sm[:,1:end .! =  m+1]
      Smj = Sm[:,Not(m+1)];
      if rank(Smj)>=n+m+1;
         XY = Smj\DT;
         Y = [XY[1:m];0];
         X = XY[m+1:end] ;
      else
         XY = Sm\DT;
         Y = XY[1:m+1];
         X = XY[m+2:end] ; 
      end
   elseif m==n-1
      XY = Sm\DT;
      Y = XY[1:m+1];
      X = XY[m+2:2*m+2] ;
   else
      ind_error = 1;
      alert("asigne_polos: Número de polos insuficiente")
   end
   K = tf(real.(X),real.(Y));
   T = feedback(K * P);
   Gur = minreal(T/P);
   S = minreal(1 - T);
   K, T, Gur, S, ind_error
end

function convol(h,u)
    n = length(h);
    m = length(u);
    l = n + m-1;
    w = zeros(Complex,l,1);
    for i in 1:l
        for j in max(1,i+1-m):min(i,n)
            w[i] = w[i] + h[j]*u[i-j+1];
        end
    end
    w;
end

function poly(raices)
   n = length(raices)
   pol = [1.0+0.0im, -raices[1]]
   for i in 2:n
       pol = convol(pol, [1.0+0.0im, -raices[i]])
   end
   return real.(pol)
end
##2
function reverse_poly(polyn)
    n = length(polyn)
    poly_r = zeros(Complex,n)
    for i in 1:n
        poly_r[i]  = polyn[n + 1 - i]
    end
    return poly_r
end


function polroots(polin)
   raices_p = PolynomialRoots.roots(reverse_poly(polin))
end

function lq_1(P,q)
   n, d = tfdata(P);
   Qd = convol(d, cambia_signo(d));
   Qn = q * convol(n, cambia_signo(n));
   Qa = zeros(Complex,length(Qd));
   l = length(Qd) - length(Qn);
   Qa[l+1:end] = Qn;
   Q = Qd .+ Qa;
   r = polroots(Q[:,1]);
   inestables = real.(r).> 0;
   indices = findall(iszero,inestables);
   dT = poly(r[indices,1])[:,1];
   T = tf(n,dT); 
   T = T / dcgain(T)[1];
   Gur = minreal(T / P) ;
   T, Gur
end
    

function tfdata(G)
   # Obtiene numerador y denominador de fcn de transferencia
   # para sistema SISO
    G = tf(G)
    n = num(G)[1,1]
    d = den(G)[1,1]
    n, d
end

function cambia_signo(n)
    n1 = n
    long = length(n1) - 1
    n2 = ones(long+1)
     for j in long:-2:1
       n2[j] = -1
    end
    return n2 .*n
end
   



function step_plt(sys, t_max)
   #=
   Grafica la respuesta al escalón del sistema sys
   =#
   rey = step(sys,t_max)
   si = stepinfo(rey)
   plot(si,legend=:bottomright)
end 

#=
Evalúa el máximo de la respuesta al escalón de sys
=#
function max_resp_step(sys)
   reu = step(sys)
   v_max = maximum(abs.(reu.y))
   v_max
end   

#=
Grafica respuesta al escalón de salida y de control
=#
function yu_stp_plt(T, Gur, tmax)
   re_y = step(T, tmax)
   sy = stepinfo(re_y)
   re_u = step(Gur, tmax)
   su = stepinfo(re_u)
   plot(plot(sy, label=L"y(t)", ylabel=L"y(t)", xlabel=L"t\, [s]", legend=:bottomright), 
   plot(su, label=L"u(t)", ylabel=L"u(t)",  xlabel=L"t\, [s]"))
end
