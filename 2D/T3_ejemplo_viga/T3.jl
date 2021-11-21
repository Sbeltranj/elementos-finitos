#JULIA 1.6.3



#=
-------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para TENSION PLANA usando elementos
      triangulares de 3 nodos
-------------------------------------------------------------------------------
DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramiento, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta=#

import XLSX
using Polynomials
using PyPlot
using LinearAlgebra
using Statistics
using SparseArrays
using PyCall 

include("t2ft_T3.jl") #Para las fuerzas superficiales y dibujar los
                      #diferentes gráficos.


close("all")          #cerrar ventanas


# %%Defino las constantes y variables, para hacer el código mas legible
X   = EF  = nodo      =  elemento = material = 1
x_  = NL1 = direccion      = lado = E  = Y   = 2
NL2 = y_  = desplazamiento = tix  = fpuntual = nu = 3
NL3 = tiy = rho = 4
tjx = espesor   = 5
tjy = 6
g   = 9.81

#Nombre archivo EXCEL
filename = "malla_refinada_v1.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")


# %% posición de los nodos:
# %%Se lee la posición de los nodos
T    = hcat(columns...)  

# xnod: fila=número del nodo, columna=coordenada X_=1 o Y_=2
xnod = T[:,x_:y_]   
nno  = length(xnod[:,1])

# %% definición de los grados de libertad
ngdl = 2*nno  
gdl  = [ [1:2:ngdl]' [2:2:ngdl]' ];    # grados de libertad (dos por nodo)
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
nmat = size(Ee,1)    # número de materiales

# %% Relación de cargas puntuales

columns, labels = XLSX.readtable(filename, "carga_punt")
T  = hcat(columns...)

ncp     = size(T,1)        # número de cargas puntuales
idxNODO = T[:,nodo]        
dirfp   = T[:,direccion];  # dirección carga
fp      = T[:,fpuntual];   # vector de fuerzas nodales equivalentes global

f = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f[gdl[idxNODO[i], dirfp[i]]] = fp[i]
end

# %% Se dibuja la malla de elementos finitos
NL1, NL2, NL3 = 1, 2, 3  # Definición de los EF triangulares por ndo

figure(1)

cgx = zeros(nef, 2) 
cgy = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef

   nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]

   plot(xnod[nod_ef,X], xnod[nod_ef,Y],
         color="k", linestyle="-")

    cgx[e] = (xnod[LaG[e,1],X] + xnod[LaG[e,2],X] + xnod[LaG[e,3],X])/3;
    cgy[e] = (xnod[LaG[e,1],Y] + xnod[LaG[e,2],Y] + xnod[LaG[e,3],Y])/3; 

   # Calculo la posición del centro de gravedad del triangulo
   text(cgx[e], cgy[e], "$e", fontsize=5, color=[1,0,0],
         horizontalalignment = "center", verticalalignment="center")
end
plot(xnod[:,X], xnod[:,Y], "b.")
axis("equal") 
title("Malla de elementos finitos")

# %% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
# equivalentes global

K   = spzeros(ngdl,ngdl)    # matriz de rigidez global como RALA (sparse)
N   = Array{Any}(undef,nef) # contenedor para las matrices de forma
B   = Array{Any}(undef,nef) # contenedor para las matrices de deformacion
idx = Array{Array{Int64}}(undef, nef,1) 
#De  = Array{Array{Int64}}(undef, nef,1) 

# matriz constitutiva del elemento para TENSION PLANA
De = [ Ee/(1 .-nue.^2)       Ee.*nue/(1 .-nue.^2)  0
       Ee.*nue/(1 .-nue.^2)  Ee/(1 .-nue.^2)       0
       0                     0                     Ee/(2 .*(1 .+nue)) ];

for e = 1:nef      # ciclo sobre todos los elementos finitos

    #% Cálculo de la matriz de rigidez del elemento e

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
                        0    c1      0    c2      0     c3
                        c1   b1      c2   b2      c3    b3 ];

    #% Calculo de la matriz de rigidez del EF e
    Ke = te[mat[e]].*B[e]'*De*B[e].*Ae;

    #% Calculo del vector de f.n.e. de fuerzas masicas del EF e (peso propio)
    fbe = -rhoe[mat[e]].*g*Ae*te[mat[e]]*[0; 1; 0; 1; 0; 1]/3;

    #% Ensamblo las contribuciones a las matrices globales
    idx[e] = [ gdl[LaG[e,NL1],:]; gdl[LaG[e,NL2],:]; gdl[LaG[e,NL3],:] ]
    K[idx[e],idx[e]] += Ke
    f[idx[e]]        += fbe

end

# %% Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

# %% # Calculo del vector de f.n.e. fte del EF e para fuerzas superficiales
columns, labels = XLSX.readtable(filename, "carga_distr")
T    = hcat(columns...)

idxNODO = T[:,nodo]
#nlcd = size(Ee,1)

el    = T[:,elemento];
lados = T[:, 2];      # número de lados con carga distribuída
tix   = T[:, tix];    # componentes de carga en las diferentes
tiy   = T[:, tiy];    # direcciones: tix, tiy, tjx, tjy
tjx   = T[:, tjx];
tjy   = T[:, tjy];

for i = 1:2

   local carga, lado, e
   e     = el[i,1]
   lado  = lados[i,1]
   
   #cargamos la función t2ft_T3
   carga = [tix[i] tiy[i] tjx[i] tjy[i]]
   fte   = t2ft_T3(xnod[LaG[e,: ],:], lado, carga, te[mat[e]])
   f[idx[e]] += fte

end


#se relacionan las restricciones
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

c = round.(Int, c)              # se convierte en Int64
c =  vec(c)                     # ahora de matrix a vector
d =  setdiff(1:ngdl,c);         # GDL desconocidos

# %% extraigo las submatrices y especifico las cantidades conocidas
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
a  = zeros(ngdl);  a[c] = ac;  a[d] = ad   # desplazamientos
q  = zeros(ngdl);  q[c] = qd;  q[d] = qc   # fuerzas nodales equivalentes


# %% Dibujo la malla de elementos finitos y las deformada de esta
delta  = reshape(a, 2, nno)'
escala = 20000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

figure(3)

for e = 1:nef
    nod_ef = LaG[e, [NL1, NL2, NL3, NL1]]
    plot(xnod[nod_ef, X], xnod[nod_ef, Y], color = "r")
    plot(xdef[nod_ef, X], xdef[nod_ef, Y], color = "b") 
         
end
title("Estructura deformada escalada ")
xlabel("x [m]")
ylabel("y [m]")
tight_layout()

# %% Se calcula para cada elemento las deformaciones y los esfuerzos
#separó memoria:
deform = zeros(3,nef)
esfuer = zeros(3,nef)

for e = 1:nef
   ae = a[idx[e]]                       # desplazamientos de los gdl del elemento e
   deform[:,e] = B[e]*ae                # calculo las deformaciones
   esfuer[:,e] = De[mat[e]]*deform[:,e] # calculo los esfuerzos
end

sx = esfuer[1,:];  sy = esfuer[2,:];  txy = esfuer[3,:]
ex = deform[1,:];  ey = deform[2,:];  gxy = deform[3,:]
ez = -(nue/Ee).*(sx+sy) 



s1   = (sx+sy)/2 + sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal maximo
s2   = (sx+sy)/2 - sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal minimo
tmax = (s1-s2)/2                                # esfuerzo cortante maximo
ang  = 0.5*atan.(2*txy, sx-sy)                  # ángulo de inclinacion de s1

#se carga la función para graficar "plot_def_esf_ang" desde t2ft_T3

plot_def_esf_ang(xnod, ex,  [], L"\epsilon_x(x,y)")
plot_def_esf_ang(xnod, ey,  [], L"\epsilon_y(x,y)")
plot_def_esf_ang(xnod, ez,  [], L"\epsilon_z(x,y)")
plot_def_esf_ang(xnod, gxy, [], L"\gamma_{xy}(x,y)")

## Se dibujan los resultados
plot_def_esf_ang(xnod, sx,  [], L"\sigma_x(x,y) [Pa]")
plot_def_esf_ang(xnod, sy,  [], L"\sigma_y(x,y) [Pa]")
plot_def_esf_ang(xnod, txy, [], L"\tau_{xy}(x,y) [Pa]")

## Se calculan y grafican para cada elemento los esfuerzos principales y
## sus direcciones

## s1, s2, taumax
plot_def_esf_ang(xnod, s1,   [ang],                  L"\sigma_1(x,y) [Pa]")
plot_def_esf_ang(xnod, s2,   [ang.+pi/2],            L"\sigma_2(x,y) [Pa]")
plot_def_esf_ang(xnod, tmax, [ang.+pi/4, ang.-pi/4], L"\tau_{max}(x,y) [Pa]")

## Calculo de los esfuerzos de von Mises
s3 = zeros(size(s1))   # s3 = zeros(size(s1)) de MATLAB
sv = sqrt.(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2)

plot_def_esf_ang(xnod, sv, [], L"\sigma_v(x,y) [Pa]")
title("Esfuerzos de von Mises (Pa)")

# %% imprimo los resultados de los desplazamientos (a), las fuerzas nodales
# equivalentes (f) y nodales de equilibrio (q)
#%% Reportamos los resultados a un libro EXCEL.

XLSX.openxlsx("resultados_T3.xlsx", mode="w") do xf
   sheet = xf[1]
   XLSX.rename!(sheet, "Desplazamientos")
   sheet["A1"] = ["Nodo ", "ux [m]", "uy [m]", "fx [N]", "fy [N]", "qx [N]", "qy [N]"]

   #Desplazamientos:
   a_ = reshape(a, 2, nno)'
   sheet["A2", dim=1] = collect(1:nno)
   sheet["B2", dim=1] = a_[:,1] 
   sheet["C2", dim=1] = a_[:,2] 

   f_ = reshape(f, 2, nno)'
   sheet["D2", dim=1] = f_[:,1] 
   sheet["E2", dim=1] = f_[:,2] 

   q_ = reshape(q, 2, nno)'
   sheet["F2", dim=1] = q_[:,1] 
   sheet["G2", dim=1] = q_[:,2] 

   XLSX.addsheet!(xf, "Esfuerzos_sx_sy_sxy")
   sheet = xf[2]     # EDIT: this works if there was only 1 sheet before. 
                     # If there were already 2 or more sheets: see comments below.

   sheet["A1"] = ["Nodo ", "sx [Pa]", "sy [Pa]", "txy [Pa]"]
   sheet["A2", dim=1] = collect(1:nef)
   sheet["B2", dim=1] = sx
   sheet["C2", dim=1] = sy 
   sheet["D2", dim=1] = txy

   XLSX.addsheet!(xf, "Deformaciones_ex_ey_gxy")
   sheet = xf[3]     # EDIT: this works if there was only 1 sheet before. 
                     # If there were already 2 or more sheets: see comments below.

   sheet["A1"] = ["Nodo ", "ex", "ey", "gxy"]
   sheet["A2", dim=1] = collect(1:nef)
   sheet["B2", dim=1] = ex
   sheet["C2", dim=1] = ey 
   sheet["D2", dim=1] = gxy

   XLSX.addsheet!(xf, "s1_s2_tmax_sv_theta")
   sheet = xf[4]     # EDIT: this works if there was only 1 sheet before. 
                     # If there were already 2 or more sheets: see comments below.

   sheet["A1"] = ["Nodo ", "s1[Pa]", "s2[Pa]", "tmax[Pa]", "sv[Pa]","theta[rad]"]
   sheet["A2", dim=1] = collect(1:nef)
   sheet["B2", dim=1] = s1
   sheet["C2", dim=1] = s2 
   sheet["D2", dim=1] = tmax
   sheet["E2", dim=1] = sv
   sheet["E2", dim=1] = ang

end

# se muestran los gráficos.
for i = 1:14

   display(figure(i))

end

gcf()

#%% Fin 