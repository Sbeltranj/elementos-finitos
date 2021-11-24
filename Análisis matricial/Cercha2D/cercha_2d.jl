# Julia 1.6.3

using LinearAlgebra


#clearconsole()

#%% Ejemplo 11.3 Uribe Escamilla
# Unidades en tonedadas y cm

E = 2040 #Elastic module ton/cm^2


ang = atand(300,400)


#barra   1    2     3    4    5
theta = [ang  0    -ang  0   -90]    # ángulo de inclinación
L     = [500  400   500  400 300]    # longitud barra
A     = [100  40    150  40   30]    # área barra

k = E * A ./ L  # rigidez de cada barra

# LaG: local a global: matriz que relaciona nodos locales y globales

LaG = [ 1 3  # (s lee la barra x va del nodo i al nodo j)
        1 4  # fila = barra
        3 2  # col1 = nodo global asociado a nodo local 1
        4 2  # col2 = nodo global asociado a nodo local 2
        3 4]

# gdl: grados de libertad
gdl = [ 1 2  # fila = nodo
        3 4  # col1 = gdl en direccion x
        5 6  # col2 = gdl en direccion y
        7 8]


#%% separo memoria
K = zeros(8, 8)  # Matrix global

T = Array{Array{Float64}}(undef, 5,1)
idx = Array{Array{Int64}}(undef, 5,1)

#%% ensamblo la matriz de rigidez global

for e = 1:5
        # saco los 4 gdls de la barra

        idx[e] = [ gdl[LaG[e,1],:];  gdl[LaG[e,2],:]]

        # matriz de rigidez local expresada en el sistema de coordenadas locales
        # para la barra e

        Kloc = [ k[e]  0 -k[e]  0
                 0     0  0     0
                -k[e]  0  k[e]  0
                 0     0  0     0 ]

        cs = cosd(theta[e])
        ss = sind(theta[e])  # sin y cos de la inclinación

        # matriz de transformacion de coordenadas para la barra e
        T[e] = [ cs  ss  0  0
                -ss  cs  0  0
                 0  0  cs  ss
                 0  0 -ss  cs]

        # se ensambla la matriz de rigidez local Ke en la matriz de rigidez global K

        K[idx[e],idx[e]] += T[e]'*Kloc*T[e] #.+ K[idx[e],idx[e]]

end
#%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)

c = [1; 2; 4]
d = [3; 5; 6; 7; 8]
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

#%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,c]; Kdc = K[d,c];
Kdd = K[d,d]; Kcd = K[c,d];

# desplazamientos para los gdls c = [1 2 4]
ac = zeros(3,1)
# fuerzas en los gdls d = [3 5 6 7 8]
fc = [0, 5*cosd(ang), 5*sind(ang), 0, -20] #ton

#%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac) #De acuerdo al manual implementa LU
qd = Kcc*ac + Kcd*ad;
#%% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(8,1);  q = zeros(8,1);  # separo la memoria
a[c] = ac;       q[c] = qd;
a[d] = ad;       #q[d] = qc = 0

#%% calculo las fuerzas axiales (fax) en cada barra

fax = Array{Array{Float64}}(undef, 5,1)

for e = 1:5 # para cada barra
   fax[e] = [-k[e] 0 k[e] 0]*T[e]*a[idx[e]]
end

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Desplazamientos:")
display(a)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Reacciones:")
println(qd)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Fuerza axial:")
display(fax)

#%%Fin
