include("fun.jl")
include("gauss_legendre.jl")
include("Bb_RM.jl")
include("Bs_RM.jl")

import XLSX
using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## Calculo de los desplazamientos verticales y ángulos de giro, las 
# reacciones, los momentos flectores y las fuerzas cortantes en una losa de
# Mindlin utilizando los elementos finitos de placa "QL9"


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

#= ## Se dibuja la malla de elementos finitos. 
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
plot(xnod[:,X], xnod[:,Y], "b.") =#

## Se cargan las funciones de forma junto con sus derivadas
# Se cargan las funciones de forma del elemento lagrangiano de 9 nodos 
# junto con sus derivadas con respecto a xi y a eta

include("funciones_forma_lagrangiano_9_nodos.jl")

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
NN = Array{Any}(undef,nef,n_gl_b,n_gl_b); # matrices de funciones de forma calculadas con n_gl_b puntos de integración
Bb = Array{Any}(undef,n_gl_b,27,n_gl_b,nef) #Array{Any}(undef,nef,n_gl_b,n_gl_b); # matrices de deformación generalizada de flexión
Bs = Array{Any}(undef,nef,n_gl_s,n_gl_s); # matrices de deformación generalizada de cortante

#s = Array{Any}(undef,n_gl_b,27,nef,2)
s = Array{Any}(undef,4,4,2,3)
## se ensambla la matriz de rigidez global K y el vector de fuerzas nodales
## equivalentes global f
for e = 1:nef      # ciclo sobre todos los elementos finitos
   xe = xnod[LaG[e,:],X];   
   ye = xnod[LaG[e,:],Y];    
    
   ## se calcula la matriz de rigidez de flexión Kb del elemento e 
   Kbe = zeros(3*nnoef,3*nnoef );
   det_Je_b = zeros(n_gl_b, n_gl_b); # Jacobianos con n_gl_b puntos de integración 

   for p = 1:n_gl_b
      for q = 1:n_gl_b
         xi_gl  = x_gl_b[p];
         eta_gl = x_gl_b[q];

         Bb[:, :, q, e] = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[1];
         det_Je_b[p,q]  = Bb_RM(xi_gl, eta_gl, xe, ye, dN_dxi, dN_deta)[2];

         # se arma la matriz de rigidez del elemento e
         Kbe +=  Bb[:, :, q, e]'*Dbg*Bb[:, :, q, e]*det_Je_b[p,q]*w_gl_b[p]*w_gl_b[q];
      end
   end
   
end




