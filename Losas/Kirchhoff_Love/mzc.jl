#Programa elaborado en JULIA 1.6.3
#Santiago Beltrán Jaramillo
#sbeltran@unal.edu.co

## Este código hace uso de la librería MatPlotlib (Python)
## instale dicha librería : https://matplotlib.org/stable/

## cargamos paquetes:
import XLSX
using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall

include("func_EF_MZC.jl")  #para los gráficos
close("all")          #cerrar ventanas

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## Defino las constantes y variables, para hacer el código mas legible
X   = EF  = nodo      =  elemento = material = 1
x_  = NL1 = direccion      = lado = E  = Y   = 2
NL2 = y_  = desplazamiento = fpuntual  = nu  = 3
NL3 = rho = 4
espesor   = NL4 = 5
fz  = 6
g   = 9.81

#Nombre archivo EXCEL
filename = "losa_rectangular_libro_solidos_efQ4.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")

## Se lee la posición de los nodos
T    = hcat(columns...)  
xnod = T[:,x_:y_] 
xnod = xnod .*1.0   
nno  = length(xnod[:,1]) #nno número de nodos

## definición de los grados de libertad
ngdl = 3*nno  
gdl  = [[1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]']    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,3)

## definición de elementos finitos con respecto a nodos de la hoja LaG_fz
# LaG: fila=número del elemento, columna=número del nodo local
columns, labels = XLSX.readtable(filename, "LaG_fz")
T = hcat(columns...)
LaG   = T[:,NL1:NL4]             # Definición de EFs respecto a nodos
LaG   = LaG*1         
nef   = size(LaG,1)
fz    = T[:,fz];                # relación de las cargas distribuidas
fz    = coalesce.(fz, 0.0)      # reemplazo los missing con ceros


## definición de las restricciones (hoja: restric): 
columns, labels = XLSX.readtable(filename, "restric")
T = hcat(columns...)
idxNODO = T[:,nodo]
dirdesp = T[:,direccion];      #dirección desplzamiento
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
d = setdiff(1:ngdl,c);         # GDL desconocidos

## Relación de cargas puntuales
columns, labels = XLSX.readtable(filename, "carga_punt")
T       = hcat(columns...)
idxNODO = T[:,nodo]
dirfp   = T[:,direccion]; # dirección carga
fp      = T[:,fpuntual];; # desplazamientos conocidos

##  Se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales
# equivalentes global "f"
f = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f[gdl[idxNODO[i], dirfp[i]]] = fp[i];
end

## Demás variables a utilizar 
columns, labels = XLSX.readtable(filename, "varios")
T       = hcat(columns...)
E          = 2.1e8   # módulo de elasticidad E
nu         = T[1,2]  # coeficiente de Poisson
rho        = T[2,2]  # densidad del material
g          = T[3,2]  # aceleración de la gravedad
t          = T[4,2]  # espesor de la losa

peso_propio = rho*g*t;  # peso propio por unidad de área

## Se dibuja la malla de elementos finitos. 
figure(1)
cg = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef
    nod_ef = LaG[e, [1, 2, 3, 4, 1]]
    
    plot(xnod[nod_ef,X], xnod[nod_ef,Y],
          color="k", linestyle="-")

    # Cálculo de la posición del centro de gravedad 
    cg[e,:] = [ mean(xnod[nod_ef,X]) mean(xnod[nod_ef,Y]) ]

    text(cg[e, X], cg[e, Y], "$e", fontsize=5, color=[1,0,0],
        horizontalalignment="center", verticalalignment="center")

end

title("Malla de elementos finitos")
plot(xnod[:,X], xnod[:,Y], "b.")

## se prepara memoria para los cálculos:
# ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K   = spzeros(ngdl,ngdl)                 # matriz de rigidez global como RALA (sparse)
idx = Array{Array{Int64}}(undef, nef,1) 
a_e = zeros(nef,1);  b_e = zeros(nef,1); # a y b de cada elemento (ancho y alto)

## matriz constitutiva
De = (E/(1-nu^2)) * [ 1  nu 0
                      nu 1  0
                      0  0  (1-nu)/2 ];

Dbe = (t^3/12)*De         # matriz constitutiva de flexión generalizada   
D = E*t^3/(12*(1-nu^2))   # rigidez a flexión de la placa  

## Calculo de Ke y fe
for e = 1:nef      # ciclo sobre todos los elementos finitos

    local a,b
    # Calculo de la matriz de rigidez del elemento e    
    x1 = xnod[LaG[e,1],X]
    x2 = xnod[LaG[e,2],X];   y2 = xnod[LaG[e,2],Y];
                             y3 = xnod[LaG[e,3],Y];
    
    a = (x2-x1)/2;  a_e[e] = a;
    b = (y3-y2)/2;  b_e[e] = b;
     
    # Calculo la matriz de rigidez Ke
    # Ke se calculo con el programa func_forma_MZC.m
    #https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/losas/Kirchhoff_Love/func_forma_MZC.m

    Ke = D/(a*b)*[ 
        b^2/a^2 - nu/5 + a^2/b^2 + 7/10          (2*nu)/5 + b^2/a^2 + 1/10          (2*nu)/5 + a^2/b^2 + 1/10     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10             b^2/a^2 - nu/10 + 1/10      a^2/(2*b^2) - (2*nu)/5 - 1/10 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5         nu/10 + b^2/(2*a^2) - 1/10         nu/10 + a^2/(2*b^2) - 1/10     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10      b^2/(2*a^2) - (2*nu)/5 - 1/10             a^2/b^2 - nu/10 + 1/10
              (2*nu)/5 + b^2/a^2 + 1/10 (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15                                 nu                  nu/10 - b^2/a^2 - 1/10     nu/15 + (2*b^2)/(3*a^2) - 1/15                                  0              1/10 - b^2/(2*a^2) - nu/10         b^2/(3*a^2) - nu/15 + 1/15                                  0           b^2/(2*a^2) - (2*nu)/5 - 1/10 (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15                                  0
              (2*nu)/5 + a^2/b^2 + 1/10                                 nu (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15           a^2/(2*b^2) - (2*nu)/5 - 1/10                                  0 (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15              1/10 - a^2/(2*b^2) - nu/10                                  0         a^2/(3*b^2) - nu/15 + 1/15                  nu/10 - a^2/b^2 - 1/10                                  0     nu/15 + (2*a^2)/(3*b^2) - 1/15
    nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10             nu/10 - b^2/a^2 - 1/10      a^2/(2*b^2) - (2*nu)/5 - 1/10         b^2/a^2 - nu/5 + a^2/b^2 + 7/10        -(2*nu)/5 - b^2/a^2 - 1/10          (2*nu)/5 + a^2/b^2 + 1/10      nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10     (2*nu)/5 - b^2/(2*a^2) + 1/10             a^2/b^2 - nu/10 + 1/10   7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5         1/10 - b^2/(2*a^2) - nu/10        nu/10 + a^2/(2*b^2) - 1/10
                 b^2/a^2 - nu/10 + 1/10     nu/15 + (2*b^2)/(3*a^2) - 1/15                                  0              -(2*nu)/5 - b^2/a^2 - 1/10 (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15                                -nu           (2*nu)/5 - b^2/(2*a^2) + 1/10 (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15                                  0              nu/10 + b^2/(2*a^2) - 1/10         b^2/(3*a^2) - nu/15 + 1/15                                  0
          a^2/(2*b^2) - (2*nu)/5 - 1/10                                  0 (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15               (2*nu)/5 + a^2/b^2 + 1/10                                -nu (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15                  nu/10 - a^2/b^2 - 1/10                                  0     nu/15 + (2*a^2)/(3*b^2) - 1/15              1/10 - a^2/(2*b^2) - nu/10                                  0         a^2/(3*b^2) - nu/15 + 1/15
7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5         1/10 - b^2/(2*a^2) - nu/10         1/10 - a^2/(2*b^2) - nu/10     nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10      (2*nu)/5 - b^2/(2*a^2) + 1/10             nu/10 - a^2/b^2 - 1/10         b^2/a^2 - nu/5 + a^2/b^2 + 7/10        -(2*nu)/5 - b^2/a^2 - 1/10        -(2*nu)/5 - a^2/b^2 - 1/10       nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10             nu/10 - b^2/a^2 - 1/10      (2*nu)/5 - a^2/(2*b^2) + 1/10
             nu/10 + b^2/(2*a^2) - 1/10         b^2/(3*a^2) - nu/15 + 1/15                                  0           (2*nu)/5 - b^2/(2*a^2) + 1/10 (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15                                  0              -(2*nu)/5 - b^2/a^2 - 1/10 (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15                                 nu                  b^2/a^2 - nu/10 + 1/10     nu/15 + (2*b^2)/(3*a^2) - 1/15                                  0
             nu/10 + a^2/(2*b^2) - 1/10                                  0         a^2/(3*b^2) - nu/15 + 1/15                  a^2/b^2 - nu/10 + 1/10                                  0     nu/15 + (2*a^2)/(3*b^2) - 1/15              -(2*nu)/5 - a^2/b^2 - 1/10                                 nu (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15           (2*nu)/5 - a^2/(2*b^2) + 1/10                                  0 (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15
    nu/5 + b^2/(2*a^2) - a^2/b^2 - 7/10      b^2/(2*a^2) - (2*nu)/5 - 1/10             nu/10 - a^2/b^2 - 1/10 7/10 - b^2/(2*a^2) - a^2/(2*b^2) - nu/5         nu/10 + b^2/(2*a^2) - 1/10         1/10 - a^2/(2*b^2) - nu/10     nu/5 - b^2/a^2 + a^2/(2*b^2) - 7/10             b^2/a^2 - nu/10 + 1/10      (2*nu)/5 - a^2/(2*b^2) + 1/10         b^2/a^2 - nu/5 + a^2/b^2 + 7/10          (2*nu)/5 + b^2/a^2 + 1/10         -(2*nu)/5 - a^2/b^2 - 1/10
          b^2/(2*a^2) - (2*nu)/5 - 1/10 (4*nu)/15 + (2*b^2)/(3*a^2) - 4/15                                  0              1/10 - b^2/(2*a^2) - nu/10         b^2/(3*a^2) - nu/15 + 1/15                                  0                  nu/10 - b^2/a^2 - 1/10     nu/15 + (2*b^2)/(3*a^2) - 1/15                                  0               (2*nu)/5 + b^2/a^2 + 1/10 (4*b^2)/(3*a^2) - (4*nu)/15 + 4/15                                -nu
                 a^2/b^2 - nu/10 + 1/10                                  0     nu/15 + (2*a^2)/(3*b^2) - 1/15              nu/10 + a^2/(2*b^2) - 1/10                                  0         a^2/(3*b^2) - nu/15 + 1/15           (2*nu)/5 - a^2/(2*b^2) + 1/10                                  0 (4*nu)/15 + (2*a^2)/(3*b^2) - 4/15              -(2*nu)/5 - a^2/b^2 - 1/10                                -nu (4*a^2)/(3*b^2) - (4*nu)/15 + 4/15 ];

      
    # Calculo del vector de fuerzas nodales equivalentes del elemento e
    # Fuerzas superficiales
    fe = 4*(fz[e] + peso_propio)*a*b*[ 1/4;  a/12;  b/12
                                       1/4; -a/12;  b/12 
                                       1/4; -a/12; -b/12
                                       1/4;  a/12; -b/12 ];
                                
    # Ensamblo las contribuciones a las matrices globales
    idx[e] = [ gdl[LaG[e,1],:]; gdl[LaG[e,2],:]; gdl[LaG[e,3],:]; gdl[LaG[e,4],:] ]
    K[idx[e],idx[e]] +=  Ke;
    f[idx[e],:]      +=  fe;
end

#%% Muestro la configuracion de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

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
aa = zeros(ngdl);  aa[c] = ac;  aa[d] = ad   # desplazamientos
q  = zeros(ngdl);  q[c]  = qd;   q[d] = qc   # fuerzas nodales equivalentes

## Dibujar deformada:
aa_ =  reshape(aa,3,nno)'
a_ = aa_[:,1]*1000
NL1, NL2, NL3, NL4 = 1,2,3,4

@pyimport matplotlib.tri as mtri
triangles = Vector{Vector{Int64}}(undef, 2*nef)

for e = 1:nef
    # se arma la matriz de correspondencia (LaG) de la malla
    triangles[2*e - 1] = LaG[e, [NL1, NL2, NL4]] .- 1
    triangles[2*e - 0] = LaG[e, [NL2, NL3, NL4]] .- 1
 end

triang = mtri.Triangulation(xnod[:,1], xnod[:,2], triangles=triangles) 

fig = figure()
esc = 0.8  #escala diagrama 
title("Estructura deformada $(esc) veces")
ax = fig.add_subplot(projection="3d")
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((2, 4, esc)) 
ax.plot_trisurf(triang, a_, cmap="jet")
#ax.set_xlabel("X[m]")
#ax.set_ylabel("Y[m]")
tight_layout() 


## Dibujo las reacciones
qq_ = reshape(q,3,nno)'
z = qq_
z_    = zeros(nno,3)

for i = 1:nno*3
    if z[i]!= 0
      z_[i] = z[i]

    elseif z[i] == 0
          z_[i]  = NaN
    end
end    

fig = figure()
ax = fig.add_subplot(1, 3, 1, projection="3d")
ax.stem(xnod[:,1], xnod[:,2], z_[:,1]*0.5, basefmt="None")
title("Reacciones Fz[kN]")

ax = fig.add_subplot(1, 3, 2, projection="3d")
ax.stem(xnod[:,1], xnod[:,2], z_[:,2]*80, basefmt="None")
title("Reacciones Mx[kN m]")

ax = fig.add_subplot(1, 3, 3, projection="3d")
ax.stem(xnod[:,1], xnod[:,2], z_[:,3]*60, basefmt="None")
title("Reacciones My[kN m]")

## Se calcula para cada elemento el vector de momentos en los puntos
## de Gauss
n_gl = 2;                          # orden de la cuadratura
x_gl = [ -sqrt(1/3); +sqrt(1/3) ]; # raices del polinomio de Legendre
sigma_b = Array{Any}(undef,nef,n_gl,n_gl)

for e = 1:nef
    local a,b
    a = a_e[e]; b = b_e[e];
   
    for i = 1:n_gl
        for j = 1:n_gl
            
            local xi, eta
            xi = x_gl[i];		eta = x_gl[j];

            # Se calcula matriz Db*B en los puntos de Gauss
            # Db_Bb se calculo con el programa func_forma_MZC.m
            Db_Bb = D/4 *[      # = Db*Bb
                (3*eta*nu*(xi - 1))/b^2 - (3*xi - 3*eta*xi)/a^2        ((3*xi - 1)*(eta - 1))/a^2                   (nu*(3*eta - 1)*(xi - 1))/b^2              (3*xi - 3*eta*xi)/a^2 - (3*eta*nu*(xi + 1))/b^2            ((3*xi + 1)*(eta - 1))/a^2             -(nu*(3*eta - 1)*(xi + 1))/b^2                   (3*xi + 3*eta*xi)/a^2 + (3*eta*nu*(xi + 1))/b^2          -((3*xi + 1)*(eta + 1))/a^2              -(nu*(3*eta + 1)*(xi + 1))/b^2                -(3*xi + 3*eta*xi)/a^2 - (3*eta*nu*(xi - 1))/b^2        -((3*xi - 1)*(eta + 1))/a^2              (nu*(3*eta + 1)*(xi - 1))/b^2
                (3*nu*xi*(eta - 1))/a^2 - (3*eta - 3*eta*xi)/b^2       (nu*(3*xi - 1)*(eta - 1))/a^2                ((3*eta - 1)*(xi - 1))/b^2                -(3*eta + 3*eta*xi)/b^2 - (3*nu*xi*(eta - 1))/a^2           (nu*(3*xi + 1)*(eta - 1))/a^2          -((3*eta - 1)*(xi + 1))/b^2                      (3*eta + 3*eta*xi)/b^2 + (3*nu*xi*(eta + 1))/a^2         -(nu*(3*xi + 1)*(eta + 1))/a^2              -((3*eta + 1)*(xi + 1))/b^2                (3*eta - 3*eta*xi)/b^2 - (3*nu*xi*(eta + 1))/a^2         -(nu*(3*xi - 1)*(eta + 1))/a^2           ((3*eta + 1)*(xi - 1))/b^2
                -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b)            -((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b)      -((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b)   ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b)                 -((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b)  ((3*eta + 1)*(eta - 1)*(nu - 1))/(2*a*b)       -((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b)                 ((3*xi - 1)*(nu - 1)*(xi + 1))/(2*a*b) ((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b)       ((nu - 1)*(3*eta^2 + 3*xi^2 - 4))/(2*a*b)                 ((3*xi + 1)*(nu - 1)*(xi - 1))/(2*a*b) -((3*eta - 1)*(eta + 1)*(nu - 1))/(2*a*b) ];
            
                sigma_b[e,i,j] = Db_Bb*aa[idx[e]];
        end
    end
end

## Se calcula para cada elemento el vector de cortantes en el centro del EF
QxQy = Array{Array{Float64}}(undef, nef,1) ;  # cortantes
xi = 0; eta = 0;     # centro del EF (punto de Gauss)

for e = 1:nef

    local a,b
    a = a_e[e]; b = b_e[e]
    
    ## QQ se calculo con el programa func_forma_MZC.m
    QQ = [ 
          -((3*eta)/4 - 3/4)/a^3 - (3*eta)/(4*a*b^2)   -(3*eta - 3)/(4*a^3) -(3*eta - 1)/(4*a*b^2) ((3*eta)/4 - 3/4)/a^3 + (3*eta)/(4*a*b^2)  -(3*eta - 3)/(4*a^3) (3*eta - 1)/(4*a*b^2) -((3*eta)/4 + 3/4)/a^3 - (3*eta)/(4*a*b^2)  (3*eta + 3)/(4*a^3) (3*eta + 1)/(4*a*b^2) ((3*eta)/4 + 3/4)/a^3 + (3*eta)/(4*a*b^2)  (3*eta + 3)/(4*a^3) -(3*eta + 1)/(4*a*b^2)
           -((3*xi)/4 - 3/4)/b^3 - (3*xi)/(4*a^2*b)    -(3*xi - 1)/(4*a^2*b)    -(3*xi - 3)/(4*b^3)   ((3*xi)/4 + 3/4)/b^3 + (3*xi)/(4*a^2*b) -(3*xi + 1)/(4*a^2*b)    (3*xi + 3)/(4*b^3)   -((3*xi)/4 + 3/4)/b^3 - (3*xi)/(4*a^2*b) (3*xi + 1)/(4*a^2*b)    (3*xi + 3)/(4*b^3)   ((3*xi)/4 - 3/4)/b^3 + (3*xi)/(4*a^2*b) (3*xi - 1)/(4*a^2*b)    -(3*xi - 3)/(4*b^3) ];
    QxQy[e] = -D*QQ*aa[idx[e]];
end


## Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1)  # número de elementos adyacentes
Mx  = zeros(nno,1)
My  = zeros(nno,1)
Mxy = zeros(nno,1)
Qx  = zeros(nno,1)
Qy  = zeros(nno,1)

A = [ 
   3^(1/2)/2 + 1            -1/2            -1/2   1 - 3^(1/2)/2
            -1/2   1 - 3^(1/2)/2   3^(1/2)/2 + 1            -1/2
   1 - 3^(1/2)/2            -1/2            -1/2   3^(1/2)/2 + 1
            -1/2   3^(1/2)/2 + 1   1 - 3^(1/2)/2            -1/2 ]
       
for e = 1:nef           

   Mx[LaG[e,:],:] .+=  A * [sigma_b[e,1,1][1]
                            sigma_b[e,1,2][1]
                            sigma_b[e,2,1][1]
                            sigma_b[e,2,2][1] ]

   My[LaG[e,:],:] .+=  A * [sigma_b[e,1,1][2]
                            sigma_b[e,1,2][2]
                            sigma_b[e,2,1][2]
                            sigma_b[e,2,2][2] ]
                                        
   Mxy[LaG[e,:],:] .+= A * [sigma_b[e,1,1][3]
                            sigma_b[e,1,2][3]
                            sigma_b[e,2,1][3]
                            sigma_b[e,2,2][3] ]

   Qx[LaG[e,:],:] .+=  QxQy[e][1];   
   Qy[LaG[e,:],:] .+=  QxQy[e][2];
                                          
   num_elem_ady[LaG[e,:],:] .+=  1;
end 

## Alisado (promedio de los momentos y cortantes en los nodos)
Mx  =  Mx./num_elem_ady;  
My  =  My./num_elem_ady;  
Mxy =  Mxy./num_elem_ady;   
Qx  =  Qx./num_elem_ady;  
Qy  =  Qy./num_elem_ady; 

## Se convierte el array-matrix en un vector para los gráficos.
Mx = vec(Mx); My = vec(My); Mxy = vec(Mxy)
Qx = vec(Qx); Qy = vec(Qy)

## Se calculan y grafican para cada elemento los momentos principales y
## sus direcciones
Mt_max = sqrt.(((Mx-My)/2).^2 + Mxy.^2) # momento torsión máximo
Mf1_xy = (Mx+My)/2 + Mt_max             # momento flector máximo
Mf2_xy = (Mx+My)/2 - Mt_max             # momento flector mínimo
ang_  = 0.5*atan.(2*Mxy, Mx-My)         # ángulo de inclinación de Mf1_xy


## Se calculan y grafican los cortantes Qx, Qy y los Qmaximos, junto con 
## su ángulo de inclinación
Q_max = hypot.(Qx, Qy)
ang   = atan.(Qy, Qx)

## se dibujan los gráficos:
#Momentos Mx, My, Mxy  
plot_mom_Q_ang(xnod,[ Mx, My, Mxy], [], [], [],
                [L"Momento Mx(kN-m/m)", L" Momento My(kN-m/m)", L"Momento Mxy(kN-m/m)"])
#Momentos principales
plot_mom_Q_ang(xnod,[ Mf1_xy, Mf2_xy, Mt_max], [ang_], [ang_.+pi/2], [ang_.+pi/4, ang_.-pi/4],
                [L"Mf1_{xy}(kN-m/m)", L"Mf2_{xy}(kN-m/m)", L"Mt_{max}(kN-m/m)"])
#Cortantes Qx, Qy, Qmax 
plot_mom_Q_ang(xnod,[Qx, Qy, Q_max], [],[],[ang],
                [L"Q_x(kN/m)", L"Q_y(kN/m)",  L"Q_{max}(kN/m)"])


#calculos wood_armer
#= include("wood_armer.jl")
Mxast_sup, Myast_sup, Mxast_inf, Myast_inf = WoodArmer(Mx, My, Mxy)

#Diseño de wood y armer:
dibujar_wood_armer(xnod,[Mxast_sup, Myast_sup, Mxast_inf, Myast_inf],
                [L"M_x(kN/m)", L"Q_y(kN/m)",  L"Q_{max}(kN/m)", L"Q_{max}(kN/m)"]) =#


## comparación solución analítica
u = 0.5; v = 1; xi = 1.25; eta = 1.5;
qdist = -10;
err = zeros(nno,1);
MEF = zeros(nno,1);
analitica = zeros(nno,1);
ww = 1;

include("cal_w.jl")
for i = 1:nno
    MEF[i] = aa_[i,ww];
    analitica[i] = calc_w(xnod[i,X], xnod[i,Y], E, nu, t, 2, 4, qdist, u, v, xi, eta);
    err[i] = abs((MEF[i]-analitica[i])/analitica[i]);

    if err[i] == NaN
        global err
        err[i] = 0
    end
end
println("Observe que al comparar ambos métodos los errores relativos máximos son:")
println(maximum(filter(!isnan,err)))
println("Es decir son extremadamente pequeños !!")