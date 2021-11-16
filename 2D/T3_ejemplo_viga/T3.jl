import XLSX
using Polynomials
using PyPlot
using LinearAlgebra
using Statistics
using SparseArrays
using PyCall

# %%Defino las constantes y variables


X   = EF  = nodo      =  elemento = material = 1
x_  = NL1 = direccion      = lado = E  = Y   = 2
NL2 = y_  = desplazamiento = tix  = fpuntual = nu = 3
NL3 = tiy = rho = 4
tjx = espesor   = 5
tjy = 6


#Nombre archivo EXCEL
filename = "malla_refinada_v1.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")

# %%Se lee la posición de los nodos
T    = hcat(columns...)  

xnod = T[:,x_:y_]   
nno  = length(xnod[:,1])

# %% definición de los grados de libertad
ngdl = 2*nno  
gdl  = [ [1:2:ngdl]' [2:2:ngdl]' ];    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,2)

# %% definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
columns, labels = XLSX.readtable(filename, "LaG_mat")
T = hcat(columns...)

LaG   = T[:,NL1:NL3]        # Definición de EFs respecto a nodos
nef   = size(LaG,1)
# %% definición de los materiales
mat   = T[:,5]

columns, labels = XLSX.readtable(filename, "prop_mat")
T  = hcat(columns...)

Ee   = T[:,   E]     # [Pa]     módulo de elasticidad
nue  = T[:,  nu]     # [-]      coeficiente de Poisson
rhoe = T[:, rho]     # [kg/m³]  densidad
te   = T[:, espesor] # [m]      espesor
nmat = size(Ee,1)


# %% Relación de cargas puntuales

columns, labels = XLSX.readtable(filename, "carga_punt")
T  = hcat(columns...)

ncp = size(T,1)
idxNODO = T[:,nodo]
dirfp   = T[:,direccion];
fp      = T[:,fpuntual];; # desplazamientos conocidos

f = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f[gdl[idxNODO[i], dirfp[i]]] = fp[i]
end


NL1, NL2, NL3 = 1, 2, 3

figure(1)
cgx = zeros(nef, 2) 
cgy = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef

   nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]

   plot(xnod[nod_ef,X], xnod[nod_ef,Y],
         color="k", linestyle="-")

    cgx[e] = (xnod[LaG[e,1],X] + xnod[LaG[e,2],X] + xnod[LaG[e,3],X])/3;
    cgy[e] = (xnod[LaG[e,1],Y] + xnod[LaG[e,2],Y] + xnod[LaG[e,3],Y])/3; 

   # Calculo la posicion del centro de gravedad del triangulo
   #cg[e,:] = [ mean(xnod[LaG[:]]) mean(xnod[LaG[e,[1 3 5 7]],Y]) ]

   text(cgx[e], cgy[e], "$e", fontsize=5, color=[1,0,0],
         horizontalalignment = "center", verticalalignment="center")
end
plot(xnod[:,X], xnod[:,Y], "b.")
# text(xnod[:,X], xnod[:,Y], num2str((1:nno)'), fontsize=16)
axis("equal") # falta tight
title("Malla de elementos finitos")

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
# equivalentes global
g = 9.81

K   = spzeros(ngdl,ngdl)        # matriz de rigidez global como RALA (sparse)
N   = Array{Any}(undef,nef) # contenedor para las matrices de forma
B   = Array{Any}(undef,nef) # contenedor para las matrices de deformacion
idx = Array{Array{Int64}}(undef, nef,1) 
De  = Array{Array{Int64}}(undef, nef,1) 

De = [ Ee/(1 .-nue.^2)       Ee.*nue/(1 .-nue.^2)  0
       Ee.*nue/(1 .-nue.^2)  Ee/(1 .-nue.^2)       0
       0                     0                     Ee/(2 .*(1 .+nue)) ];

for e = 1:nef      # ciclo sobre todos los elementos finitos
 #% Calculo de la matriz de rigidez del elemento e

    local x1, x2, x3, y1, y2, y3, Ae, b1, c1, b2, c2, b3, c3, Ke, fbe

    x1 = xnod[LaG[e,1],X];              y1 = xnod[LaG[e,1],Y];
    x2 = xnod[LaG[e,2],X];              y2 = xnod[LaG[e,2],Y];
    x3 = xnod[LaG[e,3],X];              y3 = xnod[LaG[e,3],Y];

    Ae = 0.5*det([  1 x1 y1      #Área del EF e
                    1 x2 y2
                    1 x3 y3]);               
    if Ae <= 0
        error("revise las coordenadas locales del EF $e.\n")
    end

    #% Calculo de la matriz de deformaciones B del EF e
    b1 = y2-y3;        c1 = x3-x2;
    b2 = y3-y1;        c2 = x1-x3;
    b3 = y1-y2;        c3 = x2-x1;

    B[e] = (1/(2*Ae))*[ b1    0      b2    0      b3    0 
                        0    c1      0    c2      0    c3
                        c1   b1      c2   b2      c3   b3 ];

    #% Calculo de la matriz de rigidez del EF e
    Ke = te[mat[e]].*B[e]'*De*B[e].*Ae;

    #% Calculo del vector de f.n.e. de fuerzas masicas del EF e (peso propio)
    fbe = -rhoe[mat[e]].*g*Ae*te[mat[e]]*[0; 1; 0; 1; 0; 1]/3;

    #% Ensamblo las contribuciones a las matrices globales
    idx[e] = [ gdl[LaG[e,NL1],:]; gdl[LaG[e,NL2],:]; gdl[LaG[e,NL3],:] ]
    K[idx[e],idx[e]] += Ke
    f[idx[e]]        += fbe

end

figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

include("t2ft_T3.jl")
# %% Relación de cargas puntuales

columns, labels = XLSX.readtable(filename, "carga_distr")
T    = hcat(columns...)

idxNODO = T[:,nodo]
nlcd = size(Ee,1)


ft = zeros(ngdl)   # fuerzas nodales equivalentes de cargas superficiales

el    = T[:,elemento];
lados = T[:,   2];
tix   = T[:, tix];
tiy   = T[:, tiy];
tjx   = T[:, tjx];
tjy   = T[:, tjy];

for i = 1:2

   #te    = 0.1
   e     = el[i,1]
   lado  = lados[i,1]

   carga = [tix[i] tiy[i] tjx[i] tjy[i]]
   fte   = t2ft_T3(xnod[LaG[i,: ],:], lado, carga, te[mat[e]])
   ft[idx[i]] += fte

end


columns, labels = XLSX.readtable(filename, "restric")
T = hcat(columns...)

idxNODO = T[:,nodo]
dirdesp = T[:,direccion];
ac      = T[:,desplazamiento]; # desplazamientos conocidos
#ac      = ac.*0.0

# %%Grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = length(idxNODO);  
c = zeros(n_apoyos, 1);        # GDL conocidos    

for i = 1:n_apoyos
  c[i,:] .= gdl[idxNODO[i], dirdesp[i]]
end

c = round.(Int, c)             # se convierte en Int64
c =  vec(c)                     # ahora de matrix a vector
d =  setdiff(1:ngdl,c);        # GDL desconocidos


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

#qc = zeros(size(d))  # cargas de equilibrio en nodos libres ( = 0 siempre)
#qc = round.(Int, qc)
# es lo mismo que qc = zeros(size(d)); de MATLAB

## resuelvo el sistema de ecuaciones
ad = Kdd \ (fc-Kdc*ac)        # calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd     # calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl);  a[c] = ac;  a[d] = ad   # desplazamientos
q = zeros(ngdl);  q[c] = qd;  #q[d] = qc   # fuerzas nodales equivalentes


deform = zeros(3,nef)
esfuer = zeros(3,nef)

for e = 1:nef
   ae = a[idx[e]]                       # desplazamientos de los gdl del elemento e
   deform[:,e] = B[e]*ae                # calculo las deformaciones
   esfuer[:,e] = De[mat[e]]*deform[:,e] # calculo los esfuerzos
end
sx = esfuer[1,:];  sy = esfuer[2,:];  txy = esfuer[3,:]
ex = deform[2,:];  ey = deform[2,:];  gxy = deform[3,:]
#ez = -(nue/Ee)*(sx+sy)          # deformaciones ez en tensión plana
#FALTA IMPLEMENTAR COMANDOS DE GRAFICACIÓN POR MATPLOTLIB
