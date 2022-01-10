# Programa elaborado en JULIA 1.6.3
# Santiago Beltrán Jaramillo

# Visitar: https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/axisimetrico
#=
-------------------------------------------------------------------------------
NOTA: este código SOLO es apropiado para el caso AXISIMÉTRICO usando elementos
      rectangulares serendípitos de 8 nodos
-------------------------------------------------------------------------------
DEFINICIÓN DEL PROBLEMA:
Calcule los desplazamientos y las reacciones en los empotramientos, las
deformaciones y los esfuerzos de la estructura mostrada en la figura adjunta=#

# Necesario tener instalada matplotlib de python para los gráficos
# pip install matplotlib (https://matplotlib.org/stable/)

import XLSX
using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall

ENV["MPLBACKEND"]="qt5agg"
include("gausslegendre_quad.jl")
include("funcion_axisimetrico.jl")
pygui(true)
close("all")          #cerrar ventanas

# ##Defino las constantes y variables, para hacer el código mas legible
X   = EF  = nodo      =  elemento  = material = 1
x_  = NL1 = direccion       = E    = Y   = 2
NL2 = y_  = desplazamiento  = tix  = fpuntual = nu = 3
NL3 = tiy = rho = 4
tjx = espesor =NL4  = 5
tjy = NL5 = 6
NL6, NL7, NL8 = 7,8,9
g   = 9.81

#Nombre archivo EXCEL (escriba el nombre del archivo a utilizar)
filename = "boussinesq.xlsx"

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

LaG   = T[:,NL1:NL8]        # Definición de EFs respecto a nodos
nef   = size(LaG,1)

# ## definición de los materiales
mat   = T[:,10]
columns, labels = XLSX.readtable(filename, "prop_mat")
T  = hcat(columns...)

Ee   = T[:,   E]     # [Pa]     módulo de elasticidad
nue  = T[:,  nu]     # [-]      coeficiente de Poisson
rhoe = T[:, rho]     # [kg/m³]  densidad
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

f = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

## Se dibuja la malla de elementos finitos
NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = 1, 2, 3, 4, 5, 6, 7, 8  

figure(1)

cg = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef

   nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]

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
title("Malla de elementos finitos")

## Funciones de forma serendipitas del elemento rectangular de 8 nodos:
# NOTA estas funciones de forma y sus derivadas se encontraron con el
# programa c5_funciones_forma_lagrangianos_rect_2D_8_nodos.m

Nforma(xi,eta) = [
      -((eta - 1)*(xi - 1)*(eta + xi + 1))/4       # N1
      ((xi^2 - 1)*(eta - 1))/2                     # N2
      ((eta - 1)*(xi + 1)*(eta - xi + 1))/4        # N3
      -((eta^2 - 1)*(xi + 1))/2                    # N4
      ((eta + 1)*(xi + 1)*(eta + xi - 1))/4        # N5
      -((xi^2 - 1)*(eta + 1))/2                    # N6
      ((eta + 1)*(xi - 1)*(xi - eta + 1))/4        # N7
      ((eta^2 - 1)*(xi - 1))/2                 ]   # N8

## Derivadas de N con respecto a xi
dN_dxi(xi,eta) = [
      -((eta + 2*xi)*(eta - 1))/4                  # dN1_dxi
      eta*xi - xi                                  # dN2_dxi
      ((eta - 2*xi)*(eta - 1))/4                   # dN3_dxi
      1/2 - eta^2/2                                # dN4_dxi
      ((eta + 2*xi)*(eta + 1))/4                   # dN5_dxi
      -xi*(eta + 1)                                # dN6_dxi
      -((eta - 2*xi)*(eta + 1))/4                  # dN7_dxi
      eta^2/2 - 1/2                            ]   # dN8_dxi

## Derivadas de N con respecto a eta
dN_deta(xi,eta) = [
      -((2*eta + xi)*(xi - 1))/4                   # dN1_deta
      xi^2/2 - 1/2                                 # dN2_deta
      ((xi + 1)*(2*eta - xi))/4                    # dN3_deta
      -eta*(xi + 1)                                # dN4_deta
      ((2*eta + xi)*(xi + 1))/4                    # dN5_deta
      1/2 - xi^2/2                                 # dN6_deta
      -((xi - 1)*(2*eta - xi))/4                   # dN7_deta
      eta*(xi - 1)                             ]   # dN8_deta

## párametros de la cuadratura de Gauss-Legendre
# se asumira aqui el mismo orden de la cuadratura tanto en la dirección de
# xi como en la dirección de eta
n_gl = 2    # orden de la cuadratura de Gauss-Legendre


# El comando:
x_gl, w_gl  = gausslegendre_quad(n_gl)
# calcula las raices (x_gl) y los pesos (w_gl) de polinomios de Legendre
# >> [x_gl,w_gl] = gausslegendre_quad(1)
# x_gl = 0;
# w_gl = 2;
# >> [x_gl,w_gl] = gausslegendre_quad(2)
# x_gl = [  -0.577350269189626;  0.577350269189626 ];
# w_gl = [   1.000000000000000;  1.000000000000000 ];
# >> [x_gl,w_gl] = gausslegendre_quad(3)
# x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
# w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
# >> [x_gl,w_gl] = gausslegendre_quad(4)
# x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
# w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547;
# 0.347854845137453 ];

## ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global

K  = spzeros(ngdl,ngdl)        # matriz de rigidez global como RALA (sparse)
N  = Array{Any}(undef,nef,n_gl,n_gl) # contenedor para las matrices de forma
B  = Array{Any}(undef,nef,n_gl,n_gl) # contenedor para las matrices de deformacion
De = Array{Any}(undef,nmat,1)
be = Array{Any}(undef,nmat,1)


# matriz constitutiva del elemento para el caso AXISIMETRICO
for i = 1:nmat

    global De
    De = (Ee[i]/((1+nue[i]).*(1-2*nue[i])))*
                  [1 - nue[i]  nue[i]    nue[i]   0        
                  nue[i]      1-nue[i]   nue[i]   0             
                  nue[i]      nue[i]    1-nue[i]  0             
                  0           0         0         (1-2*nue[i])/2]
    global be
    be = [0; -rhoe[i]*g]  # [kgf/m³] vector de fuerzas másicas
end


idx   = Array{Array{Int64}}(undef, nef,1) 

for e = 1:nef  # ciclo sobre todos los elementos finitos
    idx[e] = [  gdl[LaG[e,1],:];  gdl[LaG[e,2],:];
                gdl[LaG[e,3],:];  gdl[LaG[e,4],:];
                gdl[LaG[e,5],:];  gdl[LaG[e,6],:];
                gdl[LaG[e,7],:];  gdl[LaG[e,8],:] ]

       # Calculo las matrices de rigidez y el vector de fuerzas nodales
   # equivalentes del elemento
   Ke = zeros(16,16)
   fe = zeros(16)
   det_Je = zeros(n_gl,n_gl) # en esta matriz se almacenaran los Jacobianos

   for p = 1:n_gl
      for q = 1:n_gl
         local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta, xe, ye


         xi_gl  = x_gl[p]
         eta_gl = x_gl[q]

         # Se evaluan las funciones de forma en los puntos de integración
         # de Gauss-Legendre
         NNforma = Nforma(xi_gl, eta_gl)

         # Se evaluan las derivadas de las funciones de forma en los puntos
         # de integración de Gauss-Legendre
         local dx_dxi, dx_deta, dy_dxi, dy_deta, Je
         ddN_dxi  = dN_dxi( xi_gl, eta_gl);     xe = xnod[LaG[e,:],X]
         ddN_deta = dN_deta(xi_gl, eta_gl);     ye = xnod[LaG[e,:],Y]

         dx_dxi  = sum(ddN_dxi  .* xe);      dy_dxi  = sum(ddN_dxi  .* ye)
         dx_deta = sum(ddN_deta .* xe);      dy_deta = sum(ddN_deta .* ye)

         # se calcula el radio del punto de Gauss
         local r
         r = sum(NNforma .* xe)

         # Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ]

         # Se calcula el determinante del Jacobiano
         det_Je[p,q] = det(Je)

         N[e,p,q] = zeros(2,2*8)
         B[e,p,q] = zeros(4,2*8)

         for i = 1:8
            # Se ensambla la matriz de funciones de forma N
            N[e,p,q][:,[2*i-1 2*i]] = [ NNforma[i]  0
                                        0           NNforma[i] ]

            # Se ensambla la matriz de deformacion del elemento B
            dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je[p,q]
            dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je[p,q]
            B[e,p,q][:,[2*i-1 2*i]] = [ dNi_dx 0          # aqui se ensambla
                                        0      dNi_dy 
                                        NNforma[i]/r 0    # y asigna la matriz
                                        dNi_dy dNi_dx ]   # B_i
         end

         # se arma la matriz de rigidez del elemento e
         Ke += B[e,p,q]'*De*B[e,p,q]*det_Je[p,q]*r*w_gl[p]*w_gl[q]

         # vector de fuerzas nodales equivalentes
         fe += N[e,p,q]'*be*det_Je[p,q]*r*w_gl[p]*w_gl[q]


      end


   end

   Ke *= 2*pi
   fe *= 2*pi
   if any(any(det_Je .<= 0))
    error("Existen elementos con det_Je negativo en el elemento $e.\n")
 end

 K[idx[e],idx[e]] += Ke
 f[idx[e],:]      += fe
end


## Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

##  Cálculo del vector  para fuerzas superficiales
columns, labels = XLSX.readtable(filename, "carga_distr")
T    = hcat(columns...)

idxNODO = T[:,nodo]
nlcd = size(idxNODO,1)

el    = T[:,elemento];
lados = T[:, 2];      # número de lados con carga distribuída
tix   = T[:, tix];    # componentes de carga en las diferentes
tiy   = T[:, tiy];    # direcciones: tix, tiy, tjx, tjy
tjx   = T[:, tjx];
tjy   = T[:, tjy];
tkx   = T[:, 7];
tky   = T[:, 8];

for i = 1:nlcd

    local carga, lado, e
    e     = el[i,1]
    lado  = lados[i,1]
    
    #cargamos la función t2ft_T3
    carga = [tix[i] tiy[i] tjx[i] tjy[i] tkx[i] tky[i]]
    fte   = t2ft_R89(xnod[LaG[e,: ],:], lado, carga)
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
ad = Kdd \ ((fc+qc)-Kdc*ac)     # cálculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd         # cálculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl);  a[c] = ac;  a[d] = ad   # desplazamientos
q  = zeros(ngdl);  q[c] = qd;  q[d] = qc   # fuerzas nodales equivalentes



## Dibujo la malla de elementos finitos y las deformada de esta
delta  = reshape(a, 2, nno)'
escala = 1000                  # factor de escalamiento de la deformada
xdef   = xnod + escala*delta    # posición de la deformada

figure(3)

for e = 1:nef
    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
    plot(xnod[nod_ef, X], xnod[nod_ef, Y], color = "r")
    plot(xdef[nod_ef, X], xdef[nod_ef, Y], color = "b") 
         
end
title("Estructura deformada escalada $(escala) veces ")
plt.gca().set_aspect("equal", adjustable="box")
xlabel("x [m]")
ylabel("y [m]")
tight_layout()

## Se calcula para cada elemento las deformaciones y los esfuerzos
def = Array{Any}(undef,nef,n_gl,n_gl)
esf = Array{Any}(undef,nef,n_gl,n_gl)

for e = 1:nef
   ae = a[idx[e]]            # desplazamientos de los gdl del elemento e

   for pp = 1:n_gl
      for qq = 1:n_gl
         def[e,pp,qq] = B[e,pp,qq]*ae     # calculo las deformaciones
         esf[e,pp,qq] = De*def[e,pp,qq]   # calculo los esfuerzos
      end
   end
end

## Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno)  # número de elementos adyacentes
sr  = zeros(nno); sz  = zeros(nno); st  = zeros(nno); trz = zeros(nno)
er  = zeros(nno); ez  = zeros(nno); et  = zeros(nno); grz = zeros(nno)


A = [
     3^(1/2)/2 + 1            -1/2            -1/2   1 - 3^(1/2)/2
   3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4
              -1/2   1 - 3^(1/2)/2   3^(1/2)/2 + 1            -1/2
   1/4 - 3^(1/2)/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 3^(1/2)/4 + 1/4
     1 - 3^(1/2)/2            -1/2            -1/2   3^(1/2)/2 + 1
   1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4
              -1/2   3^(1/2)/2 + 1   1 - 3^(1/2)/2            -1/2
   3^(1/2)/4 + 1/4 3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4 1/4 - 3^(1/2)/4 ]

# se hace la extrapolación de los esfuerzos y las deformaciones en cada elemento
# a partir de las lecturas en los puntos de Gauss

for e = 1:nef
      sr[LaG[e,:],:] .+=   A *  [ esf[e,1,1][1]
                                  esf[e,1,2][1]
                                  esf[e,2,1][1]
                                  esf[e,2,2][1] ]
   
      sz[LaG[e,:],:] .+=   A * [ esf[e,1,1][2]
                                 esf[e,1,2][2]
                                 esf[e,2,1][2]
                                 esf[e,2,2][2] ]
   
      st[LaG[e,:],:] .+=  A *  [ esf[e,1,1][3]
                                 esf[e,1,2][3]
                                 esf[e,2,1][3]
                                 esf[e,2,2][3] ]
   
      trz[LaG[e,:],:] .+=   A *[ esf[e,1,1][4]
                                 esf[e,1,2][4]
                                 esf[e,2,1][4]
                                 esf[e,2,2][4] ]

      er[LaG[e,:],:] .+=   A * [ def[e,1,1][1]
                                 def[e,1,2][1]
                                 def[e,2,1][1]
                                 def[e,2,2][1] ]
                                 
      ez[LaG[e,:],:] .+=   A * [ def[e,1,1][2]
                                 def[e,1,2][2]
                                 def[e,2,1][2]
                                 def[e,2,2][2] ]
   
      et[LaG[e,:],:] .+=   A * [ def[e,1,1][3]
                                 def[e,1,2][3]
                                 def[e,2,1][3]
                                 def[e,2,2][3] ]

      grz[LaG[e,:],:] .+=   A *[ def[e,1,1][4]
                                 def[e,1,2][4]
                                 def[e,2,1][4]
                                 def[e,2,2][4] ]
   
      num_elem_ady[LaG[e,:],:] .+= 1
   end

## Alisado (promedio de los esfuerzos en los nodos)
sr  =  sr./num_elem_ady;  er  =  er./num_elem_ady
sz  =  sz./num_elem_ady;  ez  =  ez./num_elem_ady
st  =  st./num_elem_ady;  et  =  et./num_elem_ady
trz = trz./num_elem_ady;  grz = grz./num_elem_ady

trt = 0
ttz = 0

# %% Se calculan para cada nodo los esfuerzos principales y sus direcciones
s1 = zeros(nno);  n1 = zeros((nno, 3))
s2 = zeros(nno);  n2 = zeros((nno, 3))
s3 = zeros(nno);  n3 = zeros((nno, 3))

for i = 1:nno

    local esfppales, dirppales
    esfppales, dirppales = eigen(
                              [sr[i]   trt    trz[i]  # matriz de esfuerzos
                               trt     st[i]  ttz     # de Cauchy para 
                               trz[i]  ttz    sz[i] ]) # theta = grados

    idx_esf = reverse(sortperm(esfppales))# ordene de mayor a menor
    s1[i], s2[i], s3[i] = esfppales[idx_esf]
    n1[i,:] = dirppales[:,idx_esf[1]]
    n2[i,:] = dirppales[:,idx_esf[2]]
    n3[i,:] = dirppales[:,idx_esf[3]]
end

# Esfuerzo cortante máximo
tmax = (s1-s3)/2                               # esfuerzo cortante máximo
   
## Calculo de los esfuerzos de von Mises
sv   = sqrt.(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2)

## Gráfica del post-proceso:

# deformaciones
plot_def_esf_ang(xnod, er,  [], L"\epsilon_r")
plot_def_esf_ang(xnod, ez,  [], L"\epsilon_z")
plot_def_esf_ang(xnod, et,  [], L"\epsilon_\theta")
plot_def_esf_ang(xnod, grz, [], L"\gamma_{rz} [rad]")

# esfuerzos
plot_def_esf_ang(xnod, sr,  [], L"\sigma_r  [Pa]")
plot_def_esf_ang(xnod, sz,  [], L"\sigma_z  [Pa]")
plot_def_esf_ang(xnod, st,  [], L"\sigma_\theta  [Pa]")
plot_def_esf_ang(xnod, trz, [], L"\tau_{rz} [pa]")

# esfuerzos principales con sus orientaciones

#= plot_def_esf_ang(xnod, s1,  [ang], L"\sigma_1  [Pa]")
plot_def_esf_ang(xnod, s2,  [ang.+pi/2], L"\sigma_2  [Pa]")
plot_def_esf_ang(xnod, tmax,  [ang.+pi/4, ang.-pi/4], L"\tau_{máx}  [Pa]") =#


# esfuerzos de von Mises
#ang  = 0.5*atan.(2*trz, sr-sz) # angulo de inclinacion de s1
plot_def_esf_ang(xnod, sv, [], L"\sigma_{VM} [pa]")


XLSX.openxlsx("resultados_ejemplo_Q8_axisimetrico.xlsx", mode="w") do xf
      sheet = xf[1]
      XLSX.rename!(sheet, "Desplazamientos")
      sheet["A1"] = ["Nodo ", "ur [m]", "w [m]", "fr [N]", "fz [N]", "qr [N]", "qz [N]"]
   
      #Desplazamientos:
      a_ = reshape(a, 2, nno)'
      sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = a_[:,1]; sheet["C2", dim=1] = a_[:,2] 
      #Fuerzas 
      f_ = reshape(f, 2, nno)'; 
      sheet["D2", dim=1] = f_[:,1]; sheet["E2", dim=1] = f_[:,2] 
      # reacciones
      q_ = reshape(q, 2, nno)'
      sheet["F2", dim=1] = q_[:,1]; sheet["G2", dim=1] = q_[:,2] 
      
      XLSX.addsheet!(xf, "Esfuerzos")
      sheet = xf[2]     # EDIT: this works if there was only 1 sheet before. 
                        # If there were already 2 or more sheets: see comments below.
   
      sheet["A1"] = ["Nodo ", "sr [Pa]", "sz [Pa]", "trz [Pa]"]
      sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = sr; sheet["C2", dim=1] = sz 
      sheet["D2", dim=1] = trz
   
      XLSX.addsheet!(xf, "Deformaciones")
      sheet = xf[3]     # EDIT: this works if there was only 1 sheet before. 
                        # If there were already 2 or more sheets: see comments below.
   
      sheet["A1"] = ["Nodo ", "er", "ez", "grz"]
      sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = er; sheet["C2", dim=1] = ez 
      sheet["D2", dim=1] = grz
   
      XLSX.addsheet!(xf, "s1_s2_tmax_sv_theta")
      sheet = xf[4]     # EDIT: this works if there was only 1 sheet before. 
                        # If there were already 2 or more sheets: see comments below.
   
      sheet["A1"] = ["Nodo ", "s1[Pa]", "s2[Pa]", "tmax[Pa]", "sv[Pa]"]
      sheet["A2", dim=1] = collect(1:nno); sheet["B2", dim=1] = s1; sheet["C2", dim=1] = s2 
      sheet["D2", dim=1] = tmax; sheet["E2", dim=1] = sv; #sheet["E2", dim=1] = n1;
   end
println("Se han guardado los resultados en: resultados_ejemplo_Q8_axisimetrico.xlsx ")