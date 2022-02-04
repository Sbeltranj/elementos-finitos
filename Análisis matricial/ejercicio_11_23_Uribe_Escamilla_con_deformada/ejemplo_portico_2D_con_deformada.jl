# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/repaso_matricial/portico_2d

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

using Printf, PyPlot
ENV["MPLBACKEND"]="qt5agg"

pygui(true)
close("all")          #cerrar ventanas

#Se cargan las funciones
include("dibujar_deformada_portico.jl")
include("calc_fuerzas_nodales_equivalentes.jl")
#%% Unidades en toneladas y metros
#%% constantes
NL1 = 1; NL2 = 2; MAT = 3;
X   = 1; Y   = 2; TH  = 3;

#%% Se define la estructura
xnod = [3 4   # coordenadas de cada nodo [x, y]
        7 6 
        9 0
        0 0 ]

#% LaG: local a global: matriz que relaciona nodos locales y globales
#% fila = barra 
#% col1 = nodo global asociado a nodo local 1
#% col2 = nodo global asociado a nodo local 2
#% (se lee la barra x va del nodo i al nodo j)

#         NL1   NL2   material
barra = [   1    2    1
            4    1    2
            2    3    2 ]

LaG = barra[:, [NL1, NL2]] # local a global
mat = barra[:, MAT]         #material

#        área      inercias_y       módulo de elasticidad
#        A(m^2)     I(m^4)          E(ton/m^2)
props = [.30*.35   .30*.35^3/12     190e4
         .30*.30   .30*.30^3/12     190e4 ]

A = props[:,1];   I = props[:,2];   E = props[:,3];

nno  = size(xnod,1); # número de nodos (numero de filas de xnod)
nbar = size(LaG,1);  # número de EFs (numero de filas de LaG)
ngdl = 3*nno;        # número de grados de libertad (tres por nodo)

#%% gdl: grados de libertad
#% fila = nodo
#% col1 = gdl en dirección x
#% col2 = gdl en dirección y
#% col3 = gdl en dirección ángular antihoraria
gdl  = [ [1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]' ] # nodos vs gdl
gdl  = reshape(hcat(gdl...)',4,3)

figure(1)
#
for e = 1:nbar
    plot(xnod[LaG[e,:],X], xnod[LaG[e,:],Y], color = :blue, linewidth = 0.5)
    cgx = (xnod[LaG[e,NL1],X] + xnod[LaG[e,NL2],X])/2
    cgy = (xnod[LaG[e,NL1],Y] + xnod[LaG[e,NL2],Y])/2
    local q  = string(e)
    text(cgx, cgy, "$q", color = :red)
end
   # Calculo la posición del centro de gravedad de la barra

grid()
for n = 1:nno
  local  q  = string(n)
  text(xnod[n,X], xnod[n,Y], "$q")
  
end

gcf()
#%% cargas aplicadas (gdl carga)
cargas_aplica = [ 1.5 ]
dofs_cargados = cargas_aplica[:,1][1];


#se separa memoria 
f = zeros(ngdl, 1);
setindex!(f,dofs_cargados,gdl[1,X])

#%% fuerzas distribuidas aplicadas sobre las barras en coordenadas locales
ang1 = atan(2,4);

qxloc =   [ x -> -2.8*sin(ang1)*cos(ang1)
            x -> 0
            x -> 0]


qyloc =   [ x -> -2.8*cos(ang1)^2
            x -> 0
            x -> 0]

#%% fuerzas nodales equivalentes para las diferentes barras
# (en este ejemplo las fuerzas nodales equivalentes estas siendo 
# especificadas con respecto al sistema de coordenadas globales)

fe= Array{Array{Float64}}(undef, nbar,1)
for e = 1:nbar
    x1 = xnod[LaG[e,NL1], X];  x2 = xnod[LaG[e,NL2], X]
    y1 = xnod[LaG[e,NL1], Y];  y2 = xnod[LaG[e,NL2], Y]
    L  = hypot(x2-x1, y2-y1)

    fe[e] = calc_fuerzas_nodales_equivalentes(
         A[mat[e]], E[mat[e]], I[mat[e]], x1,x2, y1,y2, qxloc[e],qyloc[e],L)
end

#%% separó memoria
K   = zeros(ngdl,ngdl);                           # matriz de rigidez global
Ke  = Array{Array{Float64}}(undef, nbar,1)        # matriz de rigidez local en coordenadas globales
T   = Array{Array{Float64}}(undef, nbar,1)        # matriz de transformación de coordenadas
idx = Array{Array{Int64}}(undef, nbar,1)          # almacena los 6 gdls de las barras

#%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nbar  # para cada barra
   # saco los 6 gdls de la barra e
   idx[e] = [ gdl[LaG[e,1],:];  gdl[LaG[e,2],:]];

   x1 = xnod[LaG[e,NL1], X];  x2 = xnod[LaG[e,NL2], X]
   y1 = xnod[LaG[e,NL1], Y];  y2 = xnod[LaG[e,NL2], Y]

   L  =  hypot(x2-x1, y2-y1)
   c  = (x2-x1)/L;   s = (y2-y1)/L;  # seno y coseno de la inclinación

   # matriz de transformación de coordenadas para la barra e
   #c = cosd(theta[e]); s = sind(theta[e]);

   T[e] = [ c  s  0  0  0  0
           -s  c  0  0  0  0
            0  0  1  0  0  0
            0  0  0  c  s  0
            0  0  0 -s  c  0
            0  0  0  0  0  1 ]

   # matriz de rigidez local expresada en el sistema de coordenadas locales
   # para la barra e
   AE = A[mat[e]]*E[mat[e]];       L2=L^2
   EI = E[mat[e]]*I[mat[e]];       L3=L^3

   Kloc = [ AE/L   0         0        -AE/L    0          0
            0     12*EI/L3   6*EI/L2   0     -12*EI/L3   6*EI/L2
            0      6*EI/L2   4*EI/L    0      -6*EI/L2   2*EI/L
           -AE/L   0         0         AE/L    0         0
            0    -12*EI/L3  -6*EI/L2   0      12*EI/L3  -6*EI/L2
            0      6*EI/L2   2*EI/L    0      -6*EI/L2   4*EI/L];

   # matriz de rigidez local en coordenadas globales
   Ke[e]             = T[e]'*Kloc*T[e];
   K[idx[e],idx[e]] += Ke[e]
   f[idx[e]]        += fe[e] # sumo a f global
end

#%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
apoyos = [
   gdl[3,X]  0
   gdl[3,Y]  0
   gdl[3,TH] 0
   gdl[4,X]  0
   gdl[4,Y]  0
   gdl[4,TH] 0
]

c = apoyos[:,1]
d = setdiff(1:ngdl, c);

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    Recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

# %% extraigo las submatrices y especifico las cantidades conocidas

Kcc = K[c,c]; Kdc = K[d,c]; fc = f[d]
Kdd = K[d,d]; Kcd = K[c,d]; fd = f[c];

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

# desplazamientos para los gdls c 
ac = apoyos[:,2]
#%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad -fd

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);  # separó la memoria
a[c] = ac;       q[c] = qd;
a[d] = ad;      #q[d] = qc = 0


#%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
#% globales
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                                                                            ")

qe_loc= Array{Array{Float64}}(undef, nbar,1)

for e = 1:nbar # para cada barra

   qe_coord_glob = Ke[e]*a[idx[e]] - fe[e];
   print("Fuerzas internas para barra $e en coord. globales   ")
   display(qe_coord_glob)
   println("                                                                            ")
   qe_loc[e] = T[e]*qe_coord_glob
   print("Fuerzas internas para barra $e en coord. locales   ")
   display(hcat(qe_loc[e]...)')

end

vect_mov = reshape(a,3,nno)' # matrix 3x4

println("Desplazamientos nodales")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

for i = 1:nno

   u = round(1000*vect_mov[i,1], digits = 3)
   v = round(1000*vect_mov[i,2], digits = 3)

   @printf("Node %d: u = %10.3f mm v = %10.3f mm theta = %10.3f rad \n",
            i, u, v, vect_mov[i,3]) 
end

println("                                                                                    ")
qq = reshape(q,3,nno)' # matrix 3x4
println("Fuerzas nodales de equilibrio (solo imprimo los diferentes de cero)")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

for i = 1:nno
   if qq[i] != 0 && qq[i,2] != 0  && qq[i,3]  != 0

      q1 = round(qq[i,1], digits = 3)
      q2 = round(qq[i,2], digits = 3)
      @printf("Node %d: qx = %10.3f ton qy = %10.3f ton mom = %10.3f ton*m \n",
               i, q1, q2, qq[i,3]) 
   else
   end

end

#%% Dibujar la estructura y su deformada

esc_def    = 50            # escalamiento de la deformada
esc_faxial = 0.2           # escalamiento del diagrama de axiales
esc_V      = 0.3           # escalamiento del diagrama de cortantes
esc_M      = 0.3           # escalamiento del diagrama de momentos

figure(2)
title("Deformada")
xlabel("x, m")
ylabel("y, m")

figure(3)
title("Fuerza axial [ton]")
xlabel("x, m")
ylabel("y, m")

figure(4)
title("Fuerza cortante [ton]")
xlabel("x, m")
ylabel("y, m")

figure(5)
title("Momento flector [ton-m]")
xlabel("x, m")
ylabel("y, m")

for e = 1:nbar
     x1 = xnod[LaG[e,NL1], X];  x2 = xnod[LaG[e,NL2], X]
     y1 = xnod[LaG[e,NL1], Y];  y2 = xnod[LaG[e,NL2], Y]
     L  =  hypot(x2-x1, y2-y1)

     dibujar_deformada_portico(E[mat[e]],A[mat[e]],I[mat[e]],L,x1,x2,y1,y2,qxloc[e],qyloc[e],
                       T[e]*a[idx[e]],qe_loc[e], esc_def, esc_faxial, esc_V, esc_M)

end
#Fin