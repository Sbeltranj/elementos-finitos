# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/2D/ejemplo_Q6

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


#=
-------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para TENSIÓN PLANA usando elementos
      rectangulares de 4 nodos con modos incompatibles
-------------------------------------------------------------------------------
DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramientos, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta=#

# Necesario tener instalada matplotlib de python para los gráficos
# pip install matplotlib (https://matplotlib.org/stable/)

import XLSX
using Polynomials, LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK

using PlotlyJS



#ENV["MPLBACKEND"]="qt5agg"
include("gausslegendre_quad.jl")
include("funcionQ6.jl")

#pygui(true)
#close("all")          #cerrar ventanas

# ##Defino las constantes y variables, para hacer el código mas legible
X   = EF  = nodo      =  elemento  = material = 1
x_  = NL1 = direccion       = E    = Y   = 2
NL2 = y_  = desplazamiento  = tix  = fpuntual = nu = 3
NL3 = tiy = rho = 4
tjx = espesor =NL4  = 5
tjy = 6
g   = 9.81

r_ = [ 1, 2, 3, 4, 5, 6, 7, 8 ]   # GDL a retener en condensación nodal
e_ = [ 9, 10, 11, 12 ]            # GDL a eliminar en condensación nodal

#Nombre archivo EXCEL (escriba el nombre del archivo a utilizar)
#filename = "malla1.xlsx"
filename = "malla1_no_estructurada.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")

## posición de los nodos:
##Se lee la posición de los nodos
T    = hcat(columns...)  
# xnod: fila=número del nodo, columna=coordenada X_=1 o Y_=2
xnod = T[:,x_:y_]   
nno  = length(xnod[:,1])

# ## definición de los grados de libertad
ngdl = 2*nno  
gdl  = [ [1:2:ngdl]' [2:2:ngdl]' ];    # grados de libertad (dos por nodo)
gdl  = reshape(hcat(gdl...)',nno,2)

# ## definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
columns, labels = XLSX.readtable(filename, "LaG_mat")
T = hcat(columns...)

LaG   = T[:,NL1:NL4]        # Definición de EFs respecto a nodos
nef   = size(LaG,1)

# ## definición de los materiales
mat   = T[:,5]
columns, labels = XLSX.readtable(filename, "prop_mat")
T  = hcat(columns...)

Ee   = T[:,   E]     # [Pa]     módulo de elasticidad
nue  = T[:,  nu]     # [-]      coeficiente de Poisson
rhoe = T[:, rho]     # [kg/m³]  densidad
te   = T[:, espesor] # [m]      espesor
nmat = size(Ee,1)    # número de materiales

## Relación de cargas puntuales
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

## Se dibuja la malla de elementos finitos
NL1, NL2, NL3, NL4 = 1, 2, 3, 4  # Definición de los EF triangulares por ndo

#= figure(1)

cg = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef

   nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]

   plot(xnod[nod_ef,X], xnod[nod_ef,Y],
         color="b")

   cg[e,:] = [ mean(xnod[nod_ef,X]) mean(xnod[nod_ef,Y]) ]

   # Calculo la posición del centro de gravedad del triangulo
   text(cg[e, X], cg[e, Y], "$e", fontsize=10, color=[1,0,0],
         horizontalalignment="center", verticalalignment="center")
end

# en todos los nodos se dibuja un marcador y se reporta su numeración
plot(xnod[:,X], xnod[:,Y], "r*")
for i = 1:nno
    text(xnod[i, X], xnod[i, Y], "$i", color = "r")
end
axis("equal") 
title("Malla de elementos finitos") =#

## Funciones de forma y sus derivadas del elemento rectangular de 4 nodos (N1 a
#   N4) y las funciones de forma de los modos incompatibles (N5 y N6)
Nforma(xi, eta) =       [ ((eta - 1)*(xi - 1))/4    # N1
                         -((eta - 1)*(xi + 1))/4    # N2
                          ((eta + 1)*(xi + 1))/4    # N3
                         -((eta + 1)*(xi - 1))/4] # N4
#                                      1 - xi**2,    # N5
#                                      1 - eta**2  ]) # N6

# derivadas de las funciones de forma con respecto a xi
dN_dxi(xi, eta) =       [  eta/4 - 1/4         # dN1_dxi
                           1/4 - eta/4         # dN2_dxi
                           eta/4 + 1/4         # dN3_dxi
                         - eta/4 - 1/4         # dN4_dxi
                                - 2*xi         # dN5_dxi
                                     0     ]   # dN6_dxi

# derivadas de N con respecto a eta
dN_deta(xi, eta) =       [  xi/4 - 1/4         # dN1_deta
                           -xi/4 - 1/4         # dN2_deta
                            xi/4 + 1/4         # dN3_deta
                            1/4 - xi/4         # dN4_deta
                                     0         # dN5_deta
                               - 2*eta ]       # dN6_deta

## include(gauss_legendre.jl)
n_gl = 2                       # orden de la cuadratura de Gauss-Legendre
x_gl, w_gl = gausslegendre_quad(n_gl)

# se inicializan la matriz de rigidez global y los espacios en memoria que
#  almacenarán las matrices de forma y de deformación
K       = spzeros(ngdl, ngdl)         # matriz de rigidez global
inv_Kee = zeros(4, 4, nef)
Ker     = zeros(4,8, nef)

## separó memoria para los cálculos
N     = Array{Any}(undef,nef,n_gl,n_gl) # matriz de forma en cada punto de GL
B     = Array{Any}(undef,nef,n_gl,n_gl) # matriz de deformaciones en cada punto de GL
idx   = Array{Array{Int64}}(undef, nef,1) 

# matriz constitutiva del elemento para TENSION PLANA
De = [ Ee/(1 .-nue.^2)       Ee.*nue/(1 .-nue.^2)  0
       Ee.*nue/(1 .-nue.^2)  Ee/(1 .-nue.^2)       0
       0                     0                     Ee/(2 .*(1 .+nue)) ];

be = [0 -rhoe*g]  # [kgf/m³] vector de fuerzas másicas

## para cada elemento finito en la malla:
for e = 1:nef
    # se calculan con el siguiente ciclo las matrices de rigidez y el vector de
    # fuerzas nodales equivalentes del elemento usando las cuadraturas de GL
    local Ke16, fe, det_Je, Ke
    Ke16 = zeros(12, 12)
    fe   = zeros(8)
    det_Je = zeros(n_gl,n_gl)     # matriz para almacenar los jacobianos
    Ke = zeros(12,12) 
    
    idx[e] = [  gdl[LaG[e,1],:];  gdl[LaG[e,2],:];
                gdl[LaG[e,3],:];  gdl[LaG[e,4],:]]
    
    for p = 1:n_gl
        for q = 1:n_gl
            # en cada punto de la cuadratura de Gauss-Legendre se evalúan las
            # funciones de forma y sus derivadas
            local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta, xe, ye
            xi_gl  = x_gl[p]
            eta_gl = x_gl[q]

            NNforma  = Nforma(xi_gl, eta_gl)
            ddN_dxi  = dN_dxi(xi_gl, eta_gl)
            ddN_deta = dN_deta(xi_gl, eta_gl)

            # se llaman las coordenadas nodales del elemento para calcular las
            # derivadas de la función de transformación
            xe = xnod[LaG[e,:],X]
            ye = xnod[LaG[e,:],Y]
            # se calcula el Jacobiano con N1 a N4

            local dx_dxi, dx_deta, dy_dxi, dy_deta, Je
            dx_dxi  = sum(ddN_dxi[1:4] .* xe);     dy_dxi  = sum(ddN_dxi[1:4]  .* ye)
            dx_deta = sum(ddN_deta[1:4] .* xe);    dy_deta = sum(ddN_deta[1:4] .* ye)

            # con ellas se ensambla la matriz Jacobiana del elemento y se
            # calcula su determinante

            # Se ensambla la matriz Jacobiana del elemento
            Je = [ dx_dxi   dy_dxi
                   dx_deta  dy_deta ]

            det_Je[p,q] = det(Je)

            # las matrices de forma y de deformación se evalúan y se ensamblan
            # en el punto de Gauss
            N[e,p,q] = zeros(2,2*4)
            B[e,p,q] = zeros(3,2*6)

            for i = 1:4
                N[e,p,q][:,[2*i-1 2*i]] =          [NNforma[i]         0    
                                                          0          NNforma[i]]
            end

            for i = 1:6
                dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je[p,q]
                dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je[p,q]

                B[e,p,q][:,[2*i-1 2*i]] =          [dNi_dx      0 
                                                        0      dNi_dy
                                                     dNi_dy   dNi_dx]
            end
            

            # se ensamblan la matriz de rigidez del elemento y el vector de
            # fuerzas nodales equivalentes del elemento
            Ke16 += B[e,p,q]'*De*B[e,p,q].*det_Je[p,q].*te.*w_gl[p].*w_gl[q]
            
            fe   += (N[e,p,q]'.*be)[:,2].*det_Je[p,q].*te.*w_gl[p].*w_gl[q]
        end


    # se condensan los GDL jerárquicos u5, v5, u6, v6
    Krr = Ke16[r_,r_];  Ker[:, :, e] = Ke16[e_,r_]
    Kre = Ke16[r_,e_];  inv_Kee[:, :, e] = inv(Ke16[e_,e_])
    Ke  = Krr - Kre * inv_Kee[:, :, e]  * Ker[:, :, e] 

    end
    # se determina si hay puntos con jacobiano negativo, en caso tal se termina
    # el programa y se reporta

    K[idx[e],idx[e]] += Ke
    f[idx[e],:]      += fe
end

## Muestro la configuración de la matriz K (K es rala)
#= figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero") =#

##  Cálculo del vector  para fuerzas superficiales
columns, labels = XLSX.readtable(filename, "carga_distr")
T    = hcat(columns...)

idxNODO = T[:,nodo]
nlcd = size(idxNODO,1)

el    = T[:,elemento];
lados = T[:, 2];      # número de lados con carga distribuida
tix   = T[:, tix];    # componentes de carga en las diferentes
tiy   = T[:, tiy];    # direcciones: tix, tiy, tjx, tjy
tjx   = T[:, tjx];
tjy   = T[:, tjy];

for i = 1:nlcd

    local carga, lado, e
    e     = el[i,1]
    lado  = lados[i,1]
    
    #cargamos la función t2ft_T3
    carga = [tix[i] tiy[i] tjx[i] tjy[i]]
    fte   = t2ft_Q4(xnod[LaG[e,: ],:], lado, carga, te)
    f[idx[e]] += fte
 
end
## se relacionan las restricciones
columns, labels = XLSX.readtable(filename, "restric")
T = hcat(columns...)

idxNODO = T[:,nodo]
dirdesp = T[:,direccion];
ac      = T[:,desplazamiento]; # desplazamientos conocidos

# ##Grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = length(idxNODO);  
c = zeros(n_apoyos, 1);        # GDL conocidos    

for i = 1:n_apoyos
  c[i,:] .= gdl[idxNODO[i], dirdesp[i]]
end

c = round.(Int, c)              # se convierte en Int64
c =  vec(c)                     # ahora de matrix a vector
d =  setdiff(1:ngdl,c);         # GDL desconocidos

## extraigo las submatrices y especifico las cantidades conocidas
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

## Dibujo la malla de elementos finitos y las deformada de esta
delta  = reshape(a, 2, nno)'
escala = 20000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

#= figure(3)

for e = 1:nef
    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
    plot(xnod[nod_ef, X], xnod[nod_ef, Y], color = "r")
    plot(xdef[nod_ef, X], xdef[nod_ef, Y], color = "b") 
         
end
title("Estructura deformada escalada ")
xlabel("x [m]")
ylabel("y [m]")
tight_layout() =#

## Se calcula para cada elemento las deformaciones y los esfuerzos
def = Array{Any}(undef,nef,n_gl,n_gl)
esf = Array{Any}(undef,nef,n_gl,n_gl)

for e = 1:nef
    local ar, ae
    ar = a[idx[e]]            # desplazamientos de los gdl del elemento e
    ae = [ ar
         -inv_Kee[:, :, e] * Ker[:, :, e] * ar]

    for pp = 1:n_gl
       for qq = 1:n_gl
          def[e,pp,qq] = B[e,pp,qq]*ae     # calculo las deformaciones
          esf[e,pp,qq] = De*def[e,pp,qq]   # calculo los esfuerzos
       end
    end
 end

## Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno)  # número de elementos adyacentes
sx  = zeros(nno); sy  = zeros(nno); sz  = zeros(nno); txy = zeros(nno)
txz = zeros(nno); tyz = zeros(nno); ex  = zeros(nno); ey  = zeros(nno)
gxy = zeros(nno)

# matriz de extrapolación
A = [
      3^(1/2)/2 + 1             -1/2             -1/2   1 - 3^(1/2)/2
                -1/2   1 - 3^(1/2)/2   3^(1/2)/2 + 1             -1/2
      1 - 3^(1/2)/2             -1/2             -1/2   3^(1/2)/2 + 1
                -1/2   3^(1/2)/2 + 1   1 - 3^(1/2)/2            -1/2]

# se hace la extrapolación de los esfuerzos y las deformaciones en cada elemento
# a partir de las lecturas en los puntos de Gauss

for e = 1:nef
    sx[LaG[e,:],:] .+=   A *  [ esf[e,1,1][1]
                                esf[e,1,2][1]
                                esf[e,2,1][1]
                                esf[e,2,2][1] ]
 
    sy[LaG[e,:],:] .+=   A * [ esf[e,1,1][2]
                               esf[e,1,2][2]
                               esf[e,2,1][2]
                               esf[e,2,2][2] ]
 
    txy[LaG[e,:],:] .+=  A * [ esf[e,1,1][3]
                               esf[e,1,2][3]
                               esf[e,2,1][3]
                               esf[e,2,2][3] ]
 
    ex[LaG[e,:],:] .+=   A * [ def[e,1,1][1]
                               def[e,1,2][1]
                               def[e,2,1][1]
                               def[e,2,2][1] ]
 
    ey[LaG[e,:],:] .+=   A * [ def[e,1,1][2]
                               def[e,1,2][2]
                               def[e,2,1][2]
                               def[e,2,2][2] ]
 
    gxy[LaG[e,:],:] .+=  A * [ def[e,1,1][3]
                               def[e,1,2][3]
                               def[e,2,1][3]
                               def[e,2,2][3] ]
 
    num_elem_ady[LaG[e,:],:] .+= 1
 end

## Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;  ex  =  ex./num_elem_ady
sy  =  sy./num_elem_ady;  ey  =  ey./num_elem_ady
txy = txy./num_elem_ady;  gxy = gxy./num_elem_ady
## Se calculan las deformacion ez en tensión plana
ez  = -(nue/Ee).*(sx+sy);
ez  = vec(ez);

## Se calculan y grafican para cada elemento los esfuerzos principales y
## sus direcciones
# NOTA: esto solo es válido para el caso de TENSIÓN PLANA).
# En caso de DEFORMACIÓN PLANA se deben calcular los valores y vectores
# propios de la matriz de tensiónes de Cauchy
#   [dirppales[e], esfppales[e]] =  eig([sx  txy 0    # matriz de esfuerzos
#                                       txy sy  0     # de Cauchy
#                                       0   0   0]);

s1   = (sx+sy)/2 + sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal maximo
s2   = (sx+sy)/2 - sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal minimo
tmax = (s1-s2)/2                               # esfuerzo cortante maximo
ang  = 0.5*atan.(2*txy, sx-sy) # angulo de inclinacion de s1

## Calculo de los esfuerzos de von Mises
s3 = zeros(size(s1))   # s3 = zeros(size(s1)) de MATLAB
sv = sqrt.(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2);


## Meshgrid, tomado de: https://stackoverflow.com/questions/44581049/utilizing-ndgrid-meshgrid-functionality-in-julia

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

delta = 0.1
xxi  = collect(-1:delta:1)
eeta = collect(-1:delta:1)
n    = length(xxi);
xi, eta = meshgrid(xxi, eeta)
xi = xi' ; eta = eta'
w_ = 1

w = zeros(21,21);
x = zeros(21,21);
y = zeros(21,21);


pt1 = plot(1)
#gr()

for e = 1:nef
    nno = size(xnod[LaG[e,:],1],1)
    delta  = reshape(a[idx[e]], 2, nno)'
    for i = 1:21
        for j = 1:21
            w[i,j] = sum(Nforma(xi[i,j], eta[i,j]) .* delta[:,w_]);
            x[i,j] = sum(Nforma(xi[i,j], eta[i,j]) .* xnod[LaG[e,:],1]);
            y[i,j] = sum(Nforma(xi[i,j], eta[i,j]) .* xnod[LaG[e,:],2]);
        end
    end
    pt1 = plot(contour(x=x, y=y, z=w))
    
end
display(pt1)

#= #se carga la función para graficar "plot_def_esf_ang" desde funcionQ6.jl
## Se dibujan los resultados
plot_def_esf_ang(xnod, sx,  [], L"\sigma_x(x,y) [Pa]")
plot_def_esf_ang(xnod, sy,  [], L"\sigma_y(x,y) [Pa]")
plot_def_esf_ang(xnod, txy, [], L"\tau_{xy}(x,y) [Pa]")

plot_def_esf_ang(xnod, ex,  [], L"\epsilon_x(x,y)")
plot_def_esf_ang(xnod, ey,  [], L"\epsilon_y(x,y)")
plot_def_esf_ang(xnod, ez,  [], L"\epsilon_z(x,y)")
plot_def_esf_ang(xnod, gxy, [], L"\gamma_{xy}(x,y)")

## Se calculan y grafican para cada elemento los esfuerzos principales y
## sus direcciones

## s1, s2, taumax
plot_def_esf_ang(xnod, s1,   [ang],                  L"\sigma_1(x,y) [Pa]")
plot_def_esf_ang(xnod, s2,   [ang.+pi/2],            L"\sigma_2(x,y) [Pa]")
plot_def_esf_ang(xnod, tmax, [ang.+pi/4, ang.-pi/4], L"\tau_{max}(x,y) [Pa]")

## Calculo de los esfuerzos de von Mises
plot_def_esf_ang(xnod, sv, [], L"\sigma_v(x,y) [Pa]");
title("Esfuerzos de von Mises (Pa)")

## Reportamos los resultados a un libro EXCEL.

XLSX.openxlsx("resultados_Q6.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Desplazamientos")
    sheet["A1"] = ["Nodo ", "ux [m]", "uy [m]", "fx [N]", "fy [N]", "qx [N]", "qy [N]"]
 
    #Desplazamientos:
    a_ = reshape(a, 2, nno)'
    sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = a_[:,1]; sheet["C2", dim=1] = a_[:,2] 
    #Fuerzas 
    f_ = reshape(f, 2, nno)'; 
    sheet["D2", dim=1] = f_[:,1]; sheet["E2", dim=1] = f_[:,2] 
    # reacciones
    q_ = reshape(q, 2, nno)'
    sheet["F2", dim=1] = q_[:,1]; sheet["G2", dim=1] = q_[:,2] 
    
    XLSX.addsheet!(xf, "Esfuerzos_sx_sy_sxy")
    sheet = xf[2]     # EDIT: this works if there was only 1 sheet before. 
                      # If there were already 2 or more sheets: see comments below.
 
    sheet["A1"] = ["Nodo ", "sx [Pa]", "sy [Pa]", "txy [Pa]"]
    sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = sx; sheet["C2", dim=1] = sy 
    sheet["D2", dim=1] = txy
 
    XLSX.addsheet!(xf, "Deformaciones_ex_ey_gxy")
    sheet = xf[3]     # EDIT: this works if there was only 1 sheet before. 
                      # If there were already 2 or more sheets: see comments below.
 
    sheet["A1"] = ["Nodo ", "ex", "ey", "gxy"]
    sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = ex; sheet["C2", dim=1] = ey 
    sheet["D2", dim=1] = gxy
 
    XLSX.addsheet!(xf, "s1_s2_tmax_sv_theta")
    sheet = xf[4]     # EDIT: this works if there was only 1 sheet before. 
                      # If there were already 2 or more sheets: see comments below.
 
    sheet["A1"] = ["Nodo ", "s1[Pa]", "s2[Pa]", "tmax[Pa]", "sv[Pa]","theta[rad]"]
    sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = s1; sheet["C2", dim=1] = s2 
    sheet["D2", dim=1] = tmax; sheet["E2", dim=1] = sv; sheet["E2", dim=1] = ang;
end

## se reportan resultados .vtu, para visualizar en  paraview

#Consulte la documentación:
#https://jipolanco.github.io/WriteVTK.jl/dev/grids/unstructured/#Unstructured-grid

cells   = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,nef) 
for e = 1:nef
 global cells[e] = MeshCell(VTKCellTypes.VTK_QUAD,  LaG[e,:] )
end

vtkfile = vtk_grid("Q6_element", xnod[:,1].*1.0,xnod[:,2].*1.0, cells) 

vtkfile["uv"]  = a 
vtkfile["s_x"] = sx;    vtkfile["s1"] = s1;      vtkfile["ex"] = ex; vtkfile["gxy"] = gxy;
vtkfile["s_y"] = sy;    vtkfile["s2"] = s2;      vtkfile["ey"] = ey;
vtkfile["t_xy"]  = txy; vtkfile["Tmax"] = tmax;  vtkfile["ez"] = ez;

vtkfile["sv"] = sv
vtkfile["n1"] = [cos.(ang)           sin.(ang)                    ]
vtkfile["n2"] = [cos.(ang .+ pi/2)           sin.(ang .+ pi/2)    ]

outfiles = vtk_save(vtkfile)
 =#