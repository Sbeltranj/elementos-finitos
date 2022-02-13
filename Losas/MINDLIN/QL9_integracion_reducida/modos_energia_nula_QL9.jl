# Programa original (MATLAB) elaborado por:
# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/losas/Mindlin/QL9_integracion_reducida

# Traduciendo a JULIA 1.7.1 por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

## Calculo de los desplazamientos verticales y ángulos de giro, las 
# reacciones, los momentos flectores y las fuerzas cortantes en una losa de
# Mindlin utilizando los elementos finitos de placa "QL9"

#Cargamos funciones:

include("Malla_MEN.jl"); include("gauss_legendre.jl"); include("Bb_RM.jl"); include("Bs_RM.jl")
include("dibujar_QL9.jl"); include("funciones_forma_lagrangiano_9_nodos.jl"); include("dibujar_EF_Q89.jl")

## Para borrar memoria: ctrl+D --> ENTER, en Julia REPL

using Polynomials, PyPlot, SparseArrays, PyCall, WriteVTK, Statistics, LinearAlgebra

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
n_gl_b = 2; # orden de la cuadratura de GL para la integración de Kb
n_gl_s = 2; # orden de la cuadratura de GL para la integración de Ks

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
   Kbe = zeros(3*nnoef,3*nnoef );
   det_Je_b = zeros(n_gl_b, n_gl_b); # Jacobianos con n_gl_b puntos de integración 

   for p = 1:n_gl_b
      for q = 1:n_gl_b

        local xi_gl, eta_gl
        #local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta
         xi_gl  = x_gl_b[p];
         eta_gl = x_gl_b[q];

         Bb[:, :, q, e] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[1];
         det_Je_b[p,q]  = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[2];

         # se arma la matriz de rigidez del elemento e
         Kbe +=  Bb[:, :, q, e]'*Dbg*Bb[:, :, q, e]*det_Je_b[p,q]*w_gl_b[p]*w_gl_b[q];
      end
   end
   ## se calcula la matrix Ks
   Kse = zeros(3*nnoef, 3*nnoef);   
   det_Je_s = zeros(n_gl_s, n_gl_s); # Jacobianos con n_gl_s puntos de integración
   
   local xi_gl, eta_gl
   for p = 1:n_gl_s
      for q = 1:n_gl_s
        #local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta
        
         xi_gl  = x_gl_s[p];
         eta_gl = x_gl_s[q];
         
         Bs[:, :, q, e]  = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta)[1];   
         det_Je_s[p,q]   = Bs_RM(xi_gl, eta_gl, xe, ye, Nforma, dN_dxi, dN_deta)[2];

         # se arma la matriz de rigidez del elemento e
         Kse +=  Bs[:, :, q, e]'*Dsg*Bs[:, :, q, e]*det_Je_s[p,q]*w_gl_s[p]*w_gl_s[q];
                 
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

f[gdl[45,ww]] = -10
## Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

## grados de libertad del desplazamiento conocidos y desconocidos
# determino los grados de libertad correspondientes a los bordes

c = [ gdl[1,ww]; gdl[1,tx]; 
      gdl[1,ty]; gdl[37,ww];
      gdl[37,tx]; gdl[37,ty]];

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

## Se dibuja el plano medio de la malla de elementos finitos y las deformaciones de esta
#= aa_ =  reshape(a,3,nno)'
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
esc = 1.1  #escala diagrama 
title("Estructura deformada $(esc) veces")
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((2, 4, esc)) 
img = ax.plot_trisurf(triang, a_, cmap="bwr")
colorbar(img, shrink=0.7) =#

esc = 1.5
fig = plt.figure(figsize=(16, 16))
ax = plt.axes(projection="3d")
title("Modo de energía nula en el QL9")

ax.set_box_aspect((2.5, 4.5, esc)) 
#colorbar(fig, shrink=0.7)
esc = 1.1
for e = 1:nef
   dibujar_EF_Q89_RM(xnod[LaG[e,:],X], xnod[LaG[e,:],Y],Nforma, a[idx[e]]*1000, t, esc, esc);
end

fig.savefig("fig1.pdf",dpi=400)