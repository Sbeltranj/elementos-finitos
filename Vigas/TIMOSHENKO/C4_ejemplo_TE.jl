#JULIA 1.7.1

# PROGRAMA ELABORADO POR: 
# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
#https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/vigas

# Traduciendo a JULIA por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


## NOTA, me falta actualizar a modules.jl de julia, para una correcta 
## implementación y comparación con C4_ejemplo_EB.jl 


## VIGA DE TIMOSHENKO:
# Con el programa c4_func_forma_timoshenko_lineal.m se calcularon:
# * Kb = la matriz de rigidez de flexión del elemento e
# * Ks = la matriz de rigidez de cortante del elemento e
# * fe = el vector de fuerzas nodales equivalentes
# * Bb = la matriz de deformaciones por flexión
# * Bs = la matriz de deformaciones por cortante
# * N  = matriz de funciones de forma

#Cargamos paquetes:

## Para borrar memoria (ctrl+D) en consola de JULIA + ENTER.

import XLSX
using SparseArrays, PyPlot, Printf

# Programa para el cálculo de vigas de Timoshenko-E.

##Defino las constantes y variables
Y   = EF             = nodo       = 1
NL1 = TH  = tipo     = direccion  = nodos = 2
NL2 = desplazamiento = fuerza_pun = k_    = 3
E_  = 4; I_ = 5;  G_ = 6; Aas     = 7


#Nombre archivo EXCEL
filename = "viga_Uribe_Escamilla_ej_5_5.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")

## Se lee la posición de los nodos
T    = hcat(columns...)  
xnod = T[:,nodos]                      # Posición de los nodos
L    = diff(xnod)                      # Longitud de cada EF
nno  = length(xnod);                   # número de nodos
nef  = nno - 1;                        # número de elementos finitos (EF)
ngdl = 2*nno;                          # número de grados de libertad
gdl  = [ [1:2:ngdl]' [2:2:ngdl]' ];    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,2)


## Se leen la matriz de conectividad (LaG), el modulo de elasticidad, las 
# propiedades del material y las cargas
columns, labels = XLSX.readtable(filename, "LaG_EI_q")
T = hcat(columns...)

idxEF = T[:,EF]
LaG   = T[:,NL1:NL2]        # Definición de EFs respecto a nodos


E     = T[ :,E_];            # módulo de elasticidad E del EF
I     = T[ :,I_];            # momento de inercia Iz del EF
G     = T[ :,G_];            # módulo de rigidez (para viga de Timoshenko)
Aast  = T[ :,Aas];           # área de cortante (para viga de Timoshenko)
fz    = T[:,8:9];            # relación de las cargas distribuidas
fz = coalesce.(fz, 0.0)      # reemplazo los missing con ceros

##Relación de los apoyos

columns, labels = XLSX.readtable(filename, "restric")
T = hcat(columns...)

idxNODO = T[:,nodo]
dirdesp = T[:,direccion];
ac      = T[:,desplazamiento]; # desplazamientos conocidos
ac      = ac.*0.0

## Grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = length(idxNODO);  
c = zeros(n_apoyos, 1);        # GDL conocidos    

for i = 1:n_apoyos
  c[i,:] .= gdl[idxNODO[i], dirdesp[i]]
end

c = round.(Int, c)             # se convierte en Int64
c = vec(c)                     # ahora de matrix a vector
d =  setdiff(1:ngdl,c);        # GDL desconocidos

## Relación de cargas puntuales
columns, labels = XLSX.readtable(filename, "carga_punt")
T       = hcat(columns...)

idxNODO = T[:,nodo]
dirfp   = T[:,direccion];
fp      = T[:,fuerza_pun];; # desplazamientos conocidos


## Se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales
# equivalentes global "f"
f_ini = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f_ini[gdl[idxNODO[i], dirfp[i]]] = fp[i];
end

### relación de los resortes
columns, labels = XLSX.readtable(filename, "resortes")
T       = hcat(columns...)

idxNODO = T[:,nodo]
tipores = T[:,tipo]; # Y=1 (vertical), TH=2 (rotacional)
kres    = T[:,  k_]; # constante del resorte

### grados de libertad del desplazamiento conocidos y desconocidos
K_ini = spzeros(ngdl,ngdl);   # matriz de rigidez global
n_resortes = length(idxNODO);

for i = 1:n_resortes
   local idx = gdl[idxNODO[i], tipores[i]];
   K_ini[idx,idx] = kres[i];
end


## ensamblo la matriz de rigidez global y el vector de fuerzas nodales 
## equivalentes global para la viga de Timoshenko
K = K_ini;
f = f_ini;
idx   = Array{Array{Int64}}(undef, nef,1) 

for e = 1:nef     # ciclo sobre todos los elementos finitos
   
   idx[e] = [gdl[LaG[e, 1], Y ]
            gdl[LaG[e, 1], TH]
            gdl[LaG[e, 2], Y ]
            gdl[LaG[e, 2], TH]]

   EI_L = E[e]*I[e]/L[e];
   GAast_L = G[e]*Aast[e]/L[e];
   Le = L[e];
   
   # Matriz de rigidez de flexión del elemento e
   Kb = EI_L * [
      0  0  0  0
      0 +1  0 -1
      0  0  0  0
      0 -1  0 +1 ];
   
   # Matriz de rigidez de cortante del elemento e   
   global p = nef; # numero de puntos de integración (1 por EF)

   Ks = GAast_L * [   # Con cuadratura de GL de orden 1
      1     Le/2    -1     Le/2
      Le/2  Le^2/4  -Le/2  Le^2/4
     -1    -Le/2     1    -Le/2
      Le/2  Le^2/4  -Le/2  Le^2/4 ];
   
   # Use esta matriz Ks en vez de la anterior si quiere ilustrar el shear
   # locking (bloqueo de la solución). En este caso calcule la viga con 
   # h = 0.01
#=
   global p = 2*nef; # numero de puntos de integración (2 por EF)   

   Ks = GAast_L * [  # Con cuadratura de GL de orden 2
     1      Le/2    -1     Le/2
     Le/2   Le^2/3  -Le/2  Le^2/6
    -1     -Le/2     1    -Le/2
     Le/2   Le^2/6  -Le/2  Le^2/3 ];
 =#
    
   local Ke = Kb + Ks;

  # vector de fuerzas nodales equivalentes (ver Kb_Ks_timoshenko_lineal.m)
   fe = [     (Le*(2*fz[e,1] + fz[e,2]))/6
                0
              (Le*(fz[e,1] + 2*fz[e,2]))/6
                0                          ]; 
   
   K[idx[e],idx[e]] += Ke;
   f[idx[e]]        += fe;
end


## Se resuelve el sistema de ecuaciones
## extraigo las submatrices y especificó las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K[c,c]; Kdc = K[d,c]; fc = f[d]
Kdd = K[d,d]; Kcd = K[c,d]; fd = f[c];

### resuelvo el sistema de ecuaciones

ad = Kdd\(fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad -fd

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);  # separó la memoria
a[c] = ac;          q[c] = qd;          
a[d] = ad;          #q[d] = qc = 0      

## Criterio para verificar si Ks|dd es singular o invertible
j = size(Kdd,1);
s = 1;  # por gxz


if j - s*p > 0 
    println("Como j-s*p > 0, la matriz Ks|dd posiblemente es singular")
else # if j - s*p <= 0 
    println(["Como j-s*p <= 0, la matriz Ks|dd es invertible. "
          "Disminuya el numero de puntos de integración de GL o "
          "incremente el numero de EFs."])
end    

## calculo de los momentos flectores
## (se calcula el momento en el centro de cada elemento finito)
# se reserva la memoria
# recuerde que en cada elemento se calculan los momentos en las raíces de 
# los polinomios de Legendre de grado dos

xmom  = zeros(1,nef); # posición donde se calcula momento flector
mom   = zeros(1,nef); # momento flector
xib   = [ 0 ];        # raíces del polinomio de Legendre de grado 1 (vect. col)

xcor  = zeros(3,nef); # posición donde se calcula fuerza cortante
cor   = zeros(3,nef); # fuerza cortante
xis   = [ -1; 0; 1 ]; # Linea recta, pero solo crealo al centro del EF

for e = 1:nef
   Le = L[e];    
    
   # lugar donde se calcula el momento flector y la fuerza cortante
   # (centro del EF)
   xmom[:,e] = Le*xib'/2 .+ (xnod[LaG[e,1]] + xnod[LaG[e,2]])/2;
   xcor[:,e] = Le*xis'/2 .+ (xnod[LaG[e,1]] + xnod[LaG[e,2]])/2;
     
   # vector de desplazamientos nodales del elemento a^{(e)}
   ae = a[idx[e]];

   # curvatura kappa y momento flector
   Bb = [0 -1 0 1]/Le;            # matriz de deformación de flexión
   kappa = Bb*ae;                 # curvatura
   mom[:,e] = E[e]*I[e]*kappa;    # momento flector
   
   # gamma_xz y fuerza cortante
   unos = ones(size(xis));
   Bs = [ -unos/Le,  (xis.-1)/2,  unos/Le,  -(xis.+1)/2 ];
   Bs = hcat(Bs...)

   gxz = Bs*ae;                   # gamma_xz  
   cor[:,e] = -Aast[e]*G[e]*gxz;  # fuerza cortante   
end

#%% se calculan los desplazamientos al interior de cada EF
nint = 10;           # número de puntos donde se interpolara dentro del EF
xi   = collect(LinRange(-1,1,nint)); # coordenadas naturales


xx = Array{Array{Float64}}(undef, nef,1) # Interpol de posiciones (geometría) en el elemento
ww = Array{Array{Float64}}(undef, nef,1) # Interpol desplazamientos en el elemento
tt = Array{Array{Float64}}(undef, nef,1) # Interpol ángulo en el elemento


# Matriz de funciones de forma de desplazamientos y giros
Nw      = [ (1 .-xi)/2       zeros(nint,1)  (1 .+xi)/2       zeros(nint,1) ];
Nt      = [ zeros(nint,1)  (1 .-xi)/2       zeros(nint,1)  (1 .+xi)/2      ];

for e = 1:nef        # ciclo sobre todas los elementos finitos
   # vector de desplazamientos nodales del elemento a^{(e)}
   ae = a[idx[e]];

   # interpola sobre la geometría (coord naturales a geométricas)
   xx[e] = L[e]*xi/2 .+ (xnod[LaG[e,1]] + xnod[LaG[e,2]])/2;
   
   # se calcula el desplazamiento al interior del elemento finito
   ww[e] = Nw*ae;
   
   # se calcula el angulo al interior del elemento finito
   tt[e] = atan.(Nt*ae);
end



## Imprimo los resultados
println("Desplazamientos nodales")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

vect_mov = reshape(a,2,nno)' 
for i = 1:nno

   w = round(1000*vect_mov[i,1], digits = 3)
   t = round(1000*vect_mov[i,2], digits = 3)

   @printf("Node %d: w = %10.3f m  t = %10.3f rad  \n",
            i, w, t) 
end

println("                                                                                    ")

qq = reshape(q,2,nno)' # matrix 3x4


println("Fuerzas nodales de equilibrio (solo imprimo los diferentes de cero)")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

for i = 1:nno

   if qq[i] != 0 
      q1 = round(qq[i,1], digits = 3)
      q2 = round(qq[i,2], digits = 3)
      @printf("Node %3d: Ry = %10.3f KN Mz = %10.3f kN-m  \n",
               i, q1, q2) 
   else
   end
end

## Gráfico de la solución analítica y la solución por el MEF
##1) gráfico los desplazamientos de la viga


figure(1)                       # Lienzo nuevo
subplot(211)  

grid("on")                      # retícula
xlabel("Eje x [m]")             # Título del eje X
ylabel("Desplazamientos [m]")   # Título del eje Y
tight_layout()
title("Solución con el MEF para el desplazamiento ")


for e = 1:nef # ciclo sobre todos los elementos finitos
    defor = plt.plot(xx[e], ww[e], color = :red, label = nothing); 
              # gráfico solución por MEF
end

#%% 2) Gráfico de los ángulos de giro
subplot(212)
grid("on")
tight_layout()
title("Solución con el MEF para el giro ", fontsize = 12)
ylabel("Giro (rad)")             # Título del eje Y
xlabel("Eje x [m]")              # Título del eje X

for e = 1:nef # ciclo sobre todos los elementos finitos
   giro =  plt.plot(xx[e], tt[e], color = :red,
   label = nothing);             # gráfico solución por MEF
end

## 3) gráfico  de los momentos
figure(2)                        # cree un nuevo lienzo
subplot(211)

grid("on")                       # retícula
xlabel("Eje x [m]")              # Título del eje X
ylabel("Momento flector (kN-m)") # Título del eje Y
tight_layout()
title("Solución con el MEF para el momento flector")
plt.plot(xmom[:], mom[:], color = :red )


#%% 4) Gráfico de la fuerza cortante
subplot(212)

grid("on")
plt.xlabel("Eje x [m]")             # Título del eje X
plt.ylabel("Fuerza cortante (kN)")  # Título del eje Y
plt.tight_layout()
title("Solución con el MEF para la fuerza cortante")

for e = 1:nef # ciclo sobre todos los elementos finitos
    plt.plot(xcor[:,e], cor[:,e], color = :red)
    plt.plot(xcor[2,e], cor[2,e], color = :red)
end

display(figure(1))
display(figure(2))

gcf()