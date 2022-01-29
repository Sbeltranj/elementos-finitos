# Programa original (MATLAB) elaborado por:
# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/losas/Mindlin

# Traduciendo a JULIA 1.7.1 por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

## Calculo de los desplazamientos verticales y ángulos de giro, las 
# reacciones, los momentos flectores y las fuerzas cortantes en una losa de
# Mindlin utilizando los elementos finitos de placa "QQQQ-L"

#Cargamos funciones:

include("Malla1.jl"); include("gauss_legendre.jl"); include("Bb_RM.jl"); include("Bs_RM.jl")
include("dibujar_QQQQ.jl"); include("funciones_forma_lagrangiano_9_nodos.jl")
include("Bs_QQQQ_L.jl")

## Para borrar memoria: ctrl+D --> ENTER, en Julia REPL

using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK

ENV["MPLBACKEND"]="qt5agg"
pygui(true)
close("all")

## defino las variables/constantes
X = 1; Y = 2;        # un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; # lectura del código
E  = 210e9;          # módulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;            # coeficiente de Poisson
t  = 0.05;           # espesor de la losa (m)
qdistr = -10000;     # carga (N/m^2)

nef   = size(LaG,1);  # número de EFs (numero de filas de LaG)
nnoef = size(LaG,2);  # número de nodos por EF
nno   = size(xnod,1); # número de nodos (numero de filas de xnod)
ngdl  = 3*nno;        # número de grados de libertad (tres por nodo)

gdl  = [[1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]']   # nodos vs grados de libertad
gdl  = reshape(hcat(gdl...)',nno,3)

## Se dibuja la malla de elementos finitos. 
figure(1)
cg = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef
    nod_ef = LaG[e, [1, 2, 3, 4, 5, 6, 7, 8, 1]]
    
    plt.plot(xnod[nod_ef,X], xnod[nod_ef,Y],
          color="k", linestyle="-")

    # Cálculo de la posición del centro de gravedad 
    cg[e,:] = [ mean(xnod[nod_ef,X]) mean(xnod[nod_ef,Y]) ]

    plt.text(cg[e, X], cg[e, Y], "$e", fontsize=5, color=[1,0,0],
        horizontalalignment="center", verticalalignment="center")

end

title("Malla de elementos finitos")
plot(xnod[:,X], xnod[:,Y], "b.")

## Se cargan las funciones de forma junto con sus derivadas
# Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
# junto con sus derivadas con respecto a xi y a eta



## parámetros de la cuadratura de Gauss-Legendre (INTEGRACIÓN SELECTIVA)
# se asumirá aquí el mismo orden de la cuadratura tanto en la dirección de
# xi como en la dirección de eta

# se utilizara integración COMPLETA
#=
n_gl_b = 3; # orden de la cuadratura de GL para la integración de Kb
n_gl_s = 3; # orden de la cuadratura de GL para la integración de Ks
=#

# se utilizara integración SELECTIVA
n_gl_b = 3; # orden de la cuadratura de GL para la integración de Kb
n_gl_s = 3; # orden de la cuadratura de GL para la integración de Ks

# se utilizara integración REDUCIDA
#=
n_gl_b = 2; # orden de la cuadratura de GL para la integración de Kb
n_gl_s = 2; # orden de la cuadratura de GL para la integración de Ks
=#

# calcula las raíces (x_gl) y los pesos (w_gl) de polinomios de Legendre
x_gl_b, w_gl_b  = gausslegendre_quad(n_gl_b);
x_gl_s, w_gl_s  = gausslegendre_quad(n_gl_s);


## matrices constitutivas del elemento
Db = (E/(1-nu^2))* [ 1    nu   0
                     nu   1    0
                     0    0    (1-nu)/2 ]; 

G = E/(2*(1+nu));  # módulo de rigidez
alpha = 5/6;       # coeficiente de distorsión transversal de la losa de RM
Ds = diagm([alpha*G, alpha*G]);
               
Dbg = (t^3/12)*Db; # matriz constitutiva generalizada de flexión
Dsg = t*Ds;        # matriz constitutiva generalizada de cortante

## se reserva la memoria RAM de diferentes variables
K   = spzeros(ngdl,ngdl); # matriz de rigidez global como RALA (sparse)
f   = zeros(ngdl,1);     # vector de fuerzas nodales equivalentes global
idx   = Array{Array{Int64}}(undef, nef,1)     # grados de libertad de cada elemento finito

# en los siguientes contenedores se almacenara la matriz respectiva para 
# cada punto de integración: 
nno_ = length(xnod[LaG[1,:],X])*3

NN = Array{Any}(undef,nef,n_gl_b,n_gl_b); # matrices de funciones de forma calculadas con n_gl_b puntos de integración
Bb = Array{Any}(undef,3,nno_,n_gl_b,nef)  # matrices de deformación generalizada de flexión
Bs = Array{Any}(undef,2,nno_,n_gl_s,nef); # matrices de deformación generalizada de cortante

## se ensambla la matriz de rigidez global K y el vector de fuerzas nodales
## equivalentes global f
for e = 1:nef      # ciclo sobre todos los elementos finitos
   
  local xe, ye
   xe = xnod[LaG[e,:],X];   
   ye = xnod[LaG[e,:],Y];    
    
   ## se calcula la matriz de rigidez de flexión Kb del elemento e 
   Kbe = zeros(3*nnoef,3*nnoef);
   Kse = zeros(3*nnoef,3*nnoef);
   det_Je_b = zeros(n_gl_b, n_gl_b); # Jacobianos con n_gl_b puntos de integración 

   for p = 1:n_gl_b
      for q = 1:n_gl_b

        local xi_gl, eta_gl
        #local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta
         xi_gl  = x_gl_b[p];
         eta_gl = x_gl_b[q];

         Bb[:, :, q, e] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[1];
         det_Je_b[p,q]  = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[2];

         local Je_pq
         Je_pq          = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[3];
         Bs[:, :, q, e] =  Bs_QQQQ_L(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta, Je_pq)
         
         # se arma la matriz de rigidez del elemento e
         Kse +=  Bs[:, :, q, e]'*Dsg*Bs[:, :, q, e]*det_Je_b[p,q]*w_gl_s[p]*w_gl_s[q]; 
         Kbe +=  Bb[:, :, q, e]'*Dbg*Bb[:, :, q, e]*det_Je_b[p,q]*w_gl_b[p]*w_gl_b[q];
      end
   end

   ## se calcula la matriz NN
   Mbe = zeros(3*nnoef, 3*nnoef); # matriz que se utiliza en el calculo de fe   
   local xi_gl, eta_gl
   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b[p];
         eta_gl = x_gl_b[q];

         # Se evalúan las funciones de forma en los puntos de integración
         # de Gauss-Legendre
         N = Nforma(xi_gl, eta_gl);
         
         # Se ensambla la matriz de funciones de forma N
         NN[e,p,q] = zeros(3,3*nnoef);
         for i = 1:nnoef            
            NN[e,p,q][:,[3*i-2 3*i-1 3*i]] = [ 
                N[i]    0           0
                0       N[i]        0
                0       0           N[i] ];
         end
   
         # matriz requerida para calcular el vector de fuerzas nodales 
         # equivalentes (se utiliza la integración completa)
         Mbe += NN[e,p,q]'*NN[e,p,q]*det_Je_b[p,q]*w_gl_b[p]*w_gl_b[q];                                              # REVISAR !!!!!!!!!!!!!!!!   
      end
   end  
   ## se calcula el vector de fuerzas nodales equivalentes del elemento e      
   xa = xnod[LaG[e,1],X];   ya = xnod[LaG[e,1],Y];
   xb = xnod[LaG[e,5],X];   yb = xnod[LaG[e,5],Y];

   if (xa >= 0.9999 && xb <= 1.601) && (ya >= 0.9999 && yb <= 2.001)
      ffe = zeros(nnoef, 3); ffe[:,ww] .= qdistr;
      ffe = reshape(ffe', 3*nnoef,1);                                                                                             # REVISAR !!!!!!!!!!!!!
   else
      ffe = zeros(3*nnoef,1);
   end  

   fe = Mbe*ffe;   
   ## se asocian los grados de libertad del elemento locales a los globales
   idx[e] = [ gdl[LaG[e,1],:];  gdl[LaG[e,2],:];  gdl[LaG[e,3],:];  
              gdl[LaG[e,4],:];  gdl[LaG[e,5],:];  gdl[LaG[e,6],:]; 
              gdl[LaG[e,7],:];  gdl[LaG[e,8],:];  gdl[LaG[e,9],:] ];

   K[idx[e],idx[e]] += Kbe + Kse
   f[idx[e],:]      += fe
end

## Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

## grados de libertad del desplazamiento conocidos y desconocidos
# determino los grados de libertad correspondientes a los bordes
lado_x0 = findall(x -> x==0, xnod[:,X]);     lado_y0 = findall(x -> x==0, xnod[:,Y]);
lado_x2 = findall(x -> x==2, xnod[:,X]);     lado_y4 = findall(x -> x==4, xnod[:,Y]);

c = [ gdl[lado_x0,ww]; gdl[lado_x0,ty]; 
      gdl[lado_x2,ww]; gdl[lado_x2,ty];
      gdl[lado_y0,ww]; gdl[lado_y0,tx];
      gdl[lado_y4,ww]; gdl[lado_y4,tx] ];

d = setdiff(1:ngdl,c);         # GDL desconocidos

ac = zeros(length(c),1); # desplazamientos conocidos en contorno

## se extraen las submatrices y especifico las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

## extraigo las submatrices y especifico las cantidades conocidas

Kcc = K[c,c]; Kcd = K[c,d]; fd = f[c]
Kdc = K[d,c]; Kdd = K[d,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

qc = zeros(size(d))  # cargas de equilibrio en nodos libres ( = 0 siempre)
qc = round.(Int, qc)

## resuelvo el sistema de ecuaciones
ad = Kdd \ ((fc + qc)-Kdc*ac)     # cálculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd         # cálculo fuerzas de equilibrio desconocidas
a = zeros(ngdl);   a[c]  = ac;   a[d] = ad   # desplazamientos
q  = zeros(ngdl);  q[c]  = qd;   q[d] = qc   # fuerzas nodales equivalentes

## Dibujar deformada:
aa_ =  reshape(a,3,nno)'
a_ = aa_[:,1]*1000
NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = 1,2,3,4,5,6,7,8

@pyimport matplotlib.tri as mtri
triangles = Vector{Vector{Int64}}(undef, 6*nef)

for e = 1:nef

    # se arma la matriz de correspondencia (LaG) de la nueva malla triangular
    triangles[6*e - 5] = LaG[e, [NL1, NL2, NL8]] .- 1
    triangles[6*e - 4] = LaG[e, [NL2, NL3, NL4]] .- 1
    triangles[6*e - 3] = LaG[e, [NL4, NL5, NL6]] .- 1
    triangles[6*e - 2] = LaG[e, [NL2, NL4, NL6]] .- 1
    triangles[6*e - 1] = LaG[e, [NL2, NL6, NL8]] .- 1
    triangles[6*e - 0] = LaG[e, [NL6, NL7, NL8]] .- 1

end


triang = mtri.Triangulation(xnod[:,X], xnod[:,Y], triangles=triangles) 

fig = figure()
esc = 0.8  #escala diagrama 
title("Estructura deformada $(esc) veces")
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((2, 4, esc)) 
img = ax.plot_trisurf(triang, a_, cmap="bwr")
colorbar(img, shrink=0.79) 


## En los puntos de integración de Gauss-Legendre calcular:
## El vector de momentos flectores y torsores (2x2)
## El vector de fuerzas cortantes (1x1 o 2x2)
n_gl_b = 2; x_gl_b, w_gl_b  = gausslegendre_quad(n_gl_b);

# Observe que n_gl_s = 1; interpola mal la fuerza cortante.
n_gl_s = 2; x_gl_s, w_gl_s  = gausslegendre_quad(n_gl_s);

## se calcula de nuevo Bb y Bs en cada punto de GL

Bb = Array{Any}(undef,3,nno_,n_gl_b,nef)# matrices de deformación generalizada de flexión
Bs = Array{Any}(undef,2,nno_,n_gl_s,nef); # matrices de deformación generalizada de cortante

## Se calculan los momentos y las fuerzas en los puntos de GL
#sigmag_b = Array{Any}(undef,3,1,n_gl_b,nef); # momentos flectores y torsores
#sigmag_s = Array{Any}(undef,2,1,n_gl_s,nef); # fuerzas cortantes

sigmag_b = Array{Any}(undef,nef,n_gl_b,n_gl_b)
sigmag_s = Array{Any}(undef,nef,n_gl_s,n_gl_s)

for e = 1:nef      # ciclo sobre todos los elementos finitos
    local xe, ye
    xe = xnod[LaG[e,:],X];
    ye = xnod[LaG[e,:],Y];
    
    ## se calcula la matrix Bb en los puntos de integración de GL para el
    # calculo de los momentos flectores y torsores
    local xi_gl, eta_gl
    for p = 1:n_gl_b
        for q = 1:n_gl_b
            
            xi_gl  = x_gl_b[p];
            eta_gl = x_gl_b[q];
            Bb[:, :, q, e] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[1];
            sigmag_b[e,p,q] = Dbg*Bb[:, :, q, e]*a[idx[e]];
        end
    end
    
    ## se calcula la matrix Bs en los puntos de integracion de GL para el
    # calculo de las fuerzas cortantes

    local xi_gl, eta_gl
    for p = 1:n_gl_s
        for q = 1:n_gl_s
            
            xi_gl  = x_gl_s[p];
            eta_gl = x_gl_s[q];
            Bs[:, :, q, e]  = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta)[1]
            sigmag_s[e,p,q]   = Dsg*Bs[:, :, q, e]*a[idx[e]];
        end
    end
end


## Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1)  # número de elementos adyacentes
Mx  = zeros(nno,1) ; Qy  = zeros(nno,1)
My  = zeros(nno,1) ; Qx  = zeros(nno,1)
Mxy = zeros(nno,1)

# matriz de extrapolación de esfuerzos para un elemento lagrangiano de 9
# nodos
A = [ 
   3^(1/2)/2 + 1            -1/2            -1/2   1 - 3^(1/2)/2
 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4
            -1/2   1 - 3^(1/2)/2   3^(1/2)/2 + 1            -1/2
 1/4 - 3^(1/2)/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 3^(1/2)/4 + 1/4
   1 - 3^(1/2)/2            -1/2            -1/2   3^(1/2)/2 + 1
 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4
            -1/2   3^(1/2)/2 + 1   1 - 3^(1/2)/2            -1/2
 3^(1/2)/4 + 1/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 1/4 - 3^(1/2)/4
             1/4             1/4             1/4             1/4 ];

for e = 1:nef

   Mx[LaG[e,:],:] .+=  A * [sigmag_b[e,1,1][1]
                            sigmag_b[e,1,2][1]
                            sigmag_b[e,2,1][1]
                            sigmag_b[e,2,2][1] ] 


   My[LaG[e,:],:] .+=  A * [sigmag_b[e,1,1][2]
                            sigmag_b[e,1,2][2]
                            sigmag_b[e,2,1][2]
                            sigmag_b[e,2,2][2] ]
                                        
   Mxy[LaG[e,:],:] .+= A * [sigmag_b[e,1,1][3]
                            sigmag_b[e,1,2][3]
                            sigmag_b[e,2,1][3]
                            sigmag_b[e,2,2][3] ]

    
    if n_gl_s == 1

        Qx[LaG[e,:],:] += sigmag_s[e][1];
        Qy[LaG[e,:],:] += sigmag_s[e][2];

    elseif n_gl_s == 2

        Qx[LaG[e,:],:]  += A *[ sigmag_s[e,1,1][1]
                                sigmag_s[e,1,2][1]
                                sigmag_s[e,2,1][1]
                                sigmag_s[e,2,2][1] ];

        Qy[LaG[e,:],:]  += A *[ sigmag_s[e,1,1][2]
                                sigmag_s[e,1,2][2]
                                sigmag_s[e,2,1][2]
                                sigmag_s[e,2,2][2] ];
    else

    end
   num_elem_ady[LaG[e,:],:] .+=  1;

end 

## Alisado (promedio de los momentos y cortantes en los nodos)
Mx  =  Mx./num_elem_ady;  
My  =  My./num_elem_ady;  
Mxy =  Mxy./num_elem_ady;   
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady; 


## Se convierte el array-matrix en un vector para los gráficos.
Mx = vec(Mx); My = vec(My); Mxy = vec(Mxy)
Qx = vec(Qx); Qy = vec(Qy)

## Se calculan y grafican para cada elemento los momentos principales y
## sus direcciones
Mt_max = sqrt.(((Mx-My)/2).^2 + Mxy.^2) # momento torsión máximo
Mf1_xy = (Mx+My)/2 + Mt_max             # momento flector máximo
Mf2_xy = (Mx+My)/2 - Mt_max             # momento flector mínimo
ang_  = 0.5*atan.(2*Mxy, Mx-My)         # ángulo de inclinación de Mf1_xy


## Se calculan y grafican los cortantes Qx, Qy y los Qmaximos, junto con 
## su ángulo de inclinación
Q_max = hypot.(Qx, Qy)
ang   = atan.(Qy, Qx)

## se dibujan los gráficos:
#Momentos Mx, My, Mxy  
plot_mom_Q_ang(xnod,[ Mx, My, Mxy], [], [], [],
                [L"Momento Mx(kN-m/m)", L" Momento My(kN-m/m)", L"Momento Mxy(kN-m/m)"])
#Momentos principales
plot_mom_Q_ang(xnod,[ Mf1_xy, Mf2_xy, Mt_max], [ang_], [ang_.+pi/2], [ang_.+pi/4, ang_.-pi/4],
                [L"Mf1_{xy}(kN-m/m)", L"Mf2_{xy}(kN-m/m)", L"Mt_{max}(kN-m/m)"])
#Cortantes Qx, Qy, Qmax 
plot_mom_Q_ang(xnod,[Qx, Qy, Q_max], [],[],[ang],
                [L"Q_x(kN/m)", L"Q_y(kN/m)",  L"Q_{max}(kN/m)"])