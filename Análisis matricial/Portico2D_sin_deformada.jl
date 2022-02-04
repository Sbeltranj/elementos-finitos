#Julia 1.6.3

using LinearAlgebra

# Ejemplo 11.23 Uribe Escamilla sin deformada

#%% Unidades en toneladas y metros

# se definen algunas constantes que hacen el código mas legible
NL1, NL2 = 1,2
X,   Y   = 1,2

#%% defino las variables
Aviga = 0.30*0.35;       Acol  = 0.30*0.30;       # m^2    área
Iviga = 0.30*0.35^3/12;  Icol  = 0.30*0.30^3/12;  # m^4    inercia_y

#% barra   1            2            3
A     = [ Aviga        Acol         Acol         ]; # áreas
I     = [ Iviga        Icol         Icol         ]; # inercias_y
long  = [ hypot(4,2)   5            hypot(2,6)   ]; # long barra (m)
theta = [ atand(2,4)  atand(4,3)    atand(-6,2)  ]; # ángulo inclinación (grados)

#% LaG: local a global: matriz que relaciona nodos locales y globales
#% (se lee la barra x va del nodo i al nodo j)

LaG = [ 1 2         # fila = barra
        4 1         # col1 = nodo global asociado a nodo local 1
        2 3 ]       # col2 = nodo global asociado a nodo local 2

#% gdl: grados de libertad
ngdl = 12          # numero de grados de libertad

gdl  = [ 4  5  6    # fila = nodo
         7  8  9    # col1 = gdl en dirección x
         10 11 12    # col2 = gdl en dirección y
         1  2  3 ]  # col3 = gdl en dirección angular antihoraria

E = 190*10000       # ton/m^2  módulo de elasticidad

nb = size(LaG,1)    #número de barras (número de filas de LaG)
nn = size(gdl,1)    #número de nodos  (número de filas de gdl)


#%% fuerzas nodales equivalentes para las diferentes barras
# (en este ejemplo las fuerzas nodales equivalentes están siendo
# especificadas con respecto al sistema de coordenadas globales)

fe = Array{Array{Float64}}(undef, nb,1)

#        fxi    fyi    mi     fxj    fyj    mj
#        ton    ton    ton-m  ton    ton    ton-m
fe[1] = [0;     -5.60;  -3.733;  0;     -5.60;   +3.733 ]; # OJO con los signos
fe[2] = [0;      0;      0;      0;       0;      0     ]; # mirar pag 613
fe[3] = [0;      0;      0;      0;       0;      0     ];

#%% separó la memoria
K   = zeros(ngdl,ngdl);                         # matriz de rigidez global
f   = zeros(ngdl,1);                            # vector de fuerzas nodales equivalentes global
Ke  = Array{Array{Float64}}(undef, nb,1)        # matriz de rigidez local en coordenadas globales
T   = Array{Array{Float64}}(undef, nb,1)        # matriz de transformación de coordenadas
idx = Array{Array{Int64}}(undef, nb,1)          # almacena los 6 gdls de las barras

#%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nb  # para cada barra
   # saco los 6 gdls de la barra e
   idx[e] = [ gdl[LaG[e,1],:];  gdl[LaG[e,2],:]];

   # matriz de transformación de coordenadas para la barra e
   c = cosd(theta[e]); s = sind(theta[e]);

   T[e] = [ c  s  0  0  0  0
           -s  c  0  0  0  0
            0  0  1  0  0  0
            0  0  0  c  s  0
            0  0  0 -s  c  0
            0  0  0  0  0  1 ]

   # matriz de rigidez local expresada en el sistema de coordenadas locales
   # para la barra e
   AE = A[e]*E;    EI = E*I[e];    L=long[e]; L2=long[e]^2; L3=long[e]^3;

   Kloc = [ AE/L   0         0        -AE/L    0          0
            0     12*EI/L3   6*EI/L2   0     -12*EI/L3   6*EI/L2
            0      6*EI/L2   4*EI/L    0      -6*EI/L2   2*EI/L
           -AE/L   0         0         AE/L    0         0
            0    -12*EI/L3  -6*EI/L2   0      12*EI/L3  -6*EI/L2
            0      6*EI/L2   2*EI/L    0      -6*EI/L2   4*EI/L];

   # matriz de rigidez local en coordenadas globales
   Ke[e] = T[e]'*Kloc*T[e];
   K[idx[e],idx[e]] += Ke[e]
   f[idx[e]]        += fe[e] # sumo a f global
end

#%% localizó la carga puntual de 1.5 ton en el gdl 4
f[4] += 1.5
#%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
c = [1; 2; 3; 10; 11; 12]
d = setdiff(1:ngdl, c) # d = [4 5 6 7 8 9];

#%% extraigo las submatrices y especificó las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0


Kcc = K[c,c]; Kdc = K[d,c]; fc = f[d]
Kdd = K[d,d]; Kcd = K[c,d]; fd = f[c];

# desplazamientos para los gdls c = [1 2 3 10 11 12]
ac = zeros(6,1)
#%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad -fd

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);  # separó la memoria
a[c] = ac;       q[c] = qd;
a[d] = ad;       #q[d] = qc = 0

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Desplazamientos:")
display(a)
#%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
#% globales
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                                                                            ")
for e = 1:nb # para cada barra
   
   qe_coord_glob = Ke[e]*a[idx[e]] - fe[e];
   print("Fuerzas internas para barra $e en coord. globales   ")
   display(qe_coord_glob)
   println("                                                                            ")
   print("Fuerzas internas para barra $e en coord. locales   ")
   display(T[e]*qe_coord_glob)
   println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
   println("                                                                            ")

end