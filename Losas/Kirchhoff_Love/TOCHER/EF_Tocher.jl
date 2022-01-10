# Cálculo de los desplazamientos en una placa utilizando la teoría de
# Kirchhoff-Love y el elemento finito de Tocher

# Por:
# Diego Andrés Alvarez Marín
# Sebastián Jaramillo Moreno

# Traduciendo a JULIA 1.7.1 por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu,co

# cargamos paquetes:

using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall
include("func_EF_TOCHER.jl")


close("all")          #cerrar ventanas

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## defino las variables/constantes
X = 1; Y = 2; Z = 3; # un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; # lectura del código

E  = 210e9;       # módulo de elasticidad del solido (Pa) = 210GPa
nu = 0.3;         # coeficiente de Poisson
t  = 0.05;        # espesor de la losa (m)
q  = -10000;      # carga (N/m^2)

## Definimos la geometría de la losa
#Mesh_1   # malla no tan refinada
include("Mesh_2.jl") # malla muy refinada

# LaG  - definición de elementos finitos con respecto a nodos
nno  = size(xnod,1); # número de nodos (numero de filas de xnod)
ngdl = 3*nno;        # número de grados de libertad (tres por nodo)
nef  = size(LaG,1);  # número de EFs (número de filas de LaG) 

gdl  = [[1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]']    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,3)


#Relación de cargas puntuales
f = zeros(ngdl,1); # vector de fuerzas nodales equivalentes global

## Se dibuja la malla de elementos finitos. 
figure(1)
cg = zeros(nef, 2) # almacena el centro de gravedad

for e = 1:nef
    nod_ef = LaG[e, [1, 2, 3, 1]]
    
    plt.plot(xnod[nod_ef,X], xnod[nod_ef,Y],
          color="k", linestyle="-")

    # Cálculo de la posición del centro de gravedad 
    cg[e,:] = [ mean(xnod[nod_ef,X]) mean(xnod[nod_ef,Y]) ]

    plt.text(cg[e, X], cg[e, Y], "$e", fontsize=5, color=[1,0,0],
        horizontalalignment="center", verticalalignment="center")

end

plt.title("Malla de elementos finitos")
plt.plot(xnod[:,X], xnod[:,Y], "b.")


## ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K = spzeros(ngdl,ngdl); # matriz de rigidez global como RALA (sparse)

## matriz constitutiva
De = (E/(1-nu^2)) * [ 1  nu 0
                      nu 1  0
                      0  0  (1-nu)/2 ];
               
Dbe = (t^3/12)*De;       # matriz constitutiva de flexión generalizada   

D = E*t^3/(12*(1-nu^2)); # rigidez a flexión de la placa   

## matriz L
L(x,y) =    [ 0 0 0 2 0 0 6*x       2*y   0 
			  0 0 0 0 0 2   0       2*x 6*y 
			  0 0 0 0 2 0   0 4*x + 4*y   0];

## vector p
p(x,y) = [ 1; x;	y; x^2; x*y; y^2; x^3; (x^2*y+x*y^2); y^3];
            
## se define el orden de la cuadratura
include("TriGaussPoints.jl")
orden = 3;
GP = TriGaussPoints(orden); nGP = size(GP,1);
L2 = GP[:,1];   L3 = GP[:,2];   Wi = GP[:,3]/2;

## Calculo de Ke y fe
inv_A = Array{Array{Float64}}(undef, nef,1) 

for e = 1:nef      # ciclo sobre todos los elementos finitos  
    
    local x1, x2, x3, y1, y2, y3, Ae
    ## Calculo de la matriz de rigidez Ke
    x1 = xnod[LaG[e,1],X];              y1 = xnod[LaG[e,1],Y];
    x2 = xnod[LaG[e,2],X];              y2 = xnod[LaG[e,2],Y];
    x3 = xnod[LaG[e,3],X];              y3 = xnod[LaG[e,3],Y];
    
    Ae = 0.5*det([  1 x1 y1      #Área del EF e
                    1 x2 y2
                    1 x3 y3]);               
    if Ae <= 0
        error("revise las coordenadas locales del EF $e.\n")
    end

    iint_LT_Db_L_dA = zeros(9,9);

    for i = 1:nGP
        x = (1-L2[i]-L3[i])*x1 + L2[i]*x2 + L3[i]*x3;
        y = (1-L2[i]-L3[i])*y1 + L2[i]*y2 + L3[i]*y3;
        Le = L(x,y);
        iint_LT_Db_L_dA +=  Le'*Dbe*Le*Wi[i];
    end


    iint_LT_Db_L_dA = 2*Ae*iint_LT_Db_L_dA;
    
    inv_A[e] = inv(
          [ 1 x1 y1 x1^2 x1*y1 y1^2   x1^3 x1^2*y1 + x1*y1^2   y1^3
            0  1  0 2*x1    y1    0 3*x1^2    y1^2 + 2*x1*y1      0
            0  0  1    0    x1 2*y1      0    x1^2 + 2*y1*x1 3*y1^2
            1 x2 y2 x2^2 x2*y2 y2^2   x2^3 x2^2*y2 + x2*y2^2   y2^3
            0  1  0 2*x2    y2    0 3*x2^2    y2^2 + 2*x2*y2      0
            0  0  1    0    x2 2*y2      0    x2^2 + 2*y2*x2 3*y2^2
            1 x3 y3 x3^2 x3*y3 y3^2   x3^3 x3^2*y3 + x3*y3^2   y3^3
            0  1  0 2*x3    y3    0 3*x3^2    y3^2 + 2*x3*y3      0
            0  0  1    0    x3 2*y3      0    x3^2 + 2*y3*x3 3*y3^2 ]);        

    # Cálculo la matriz de rigidez Ke
    Ke = inv_A[e]'*iint_LT_Db_L_dA*inv_A[e];

    ## Cálculo del vector de fuerzas nodales equivalentes fe

    if        ((x1 >= 0.9999 && x1 <= 1.501) 
            && (x2 >= 0.9999 && x2 <= 1.501) 
            && (x3 >= 0.9999 && x3 <= 1.501) 
            && (y1 >= 0.9999 && y1 <= 2.001) 
            && (y2 >= 0.9999 && y2 <= 2.001) 
            && (y3 >= 0.9999 && y3 <= 2.001))
        int_p_dA = 0;
        for i = 1:nGP
            x = (1-L2[i]-L3[i])*x1 + L2[i]*x2 + L3[i]*x3;
            y = (1-L2[i]-L3[i])*y1 + L2[i]*y2 + L3[i]*y3;
            int_p_dA = int_p_dA .+ p(x,y)*Wi[i];
        end
        fe = inv_A[e]'*q*(2*Ae*int_p_dA);
    else
        fe = zeros(9,1);
    end        
    # Ensamblo las contribuciones a las matrices globales
    idx = [ gdl[LaG[e,1],:]; gdl[LaG[e,2],:]; gdl[LaG[e,3],:] ]
    K[idx,idx] +=  Ke;
    f[idx,:]   +=  fe;


end

## Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

## grados de libertad del desplazamiento conocidos y desconocidos
# determino los grados de libertad correspondientes a los bordes
lado_x0 = findall(x -> x==0, xnod[:,X]);     lado_y0 = findall(x -> x==0, xnod[:,Y]);
lado_x2 = findall(x -> x==2, xnod[:,X]);     lado_y4 = findall(x -> x==4, xnod[:,Y]);

c = [ gdl[lado_x0,ww]; gdl[lado_x0,ty]; 
      gdl[lado_x2,ww]; gdl[lado_x2,ty];
      gdl[lado_y0,ww]; gdl[lado_y0,tx];
      gdl[lado_y4,ww]; gdl[lado_y4,tx] ];

d = setdiff(1:ngdl,c);         # GDL desconocidos

ac = zeros(length(c),1); # desplazamientos conocidos en contorno
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
vect_mov =  reshape(aa,3,nno)'
vect_mov = vect_mov[:,1]*1000
NL1, NL2, NL3, NL4 = 1,2,3,4

@pyimport matplotlib.tri as mtri
triangles = Vector{Vector{Int64}}(undef, 2*nef)

for e = 1:nef
    # se arma la matriz de correspondencia (LaG) de la malla
    triangles[2*e - 1] = LaG[e, [NL1, NL2, NL3]] .- 1
    triangles[2*e - 0] = LaG[e, [NL2, NL3, NL3]] .- 1
 end

triang = mtri.Triangulation(xnod[:,1], xnod[:,2], triangles=triangles) 

fig = figure()
esc = 0.8  #escala diagrama 
title("Estructura deformada $(esc) veces")
ax = fig.add_subplot(projection="3d")
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((2, 4, esc)) 
ax.plot_trisurf(triang, vect_mov, cmap="jet")
plt.tight_layout() 


## Se calcula para cada elemento el vector de momentos y cortantes en
## los puntos de Gauss

sigma_b = Array{Array{Float64}}(undef, nef,nGP);   # momentos flectores Mx y My y torsor Mxy
QxQy    = Array{Array{Float64}}(undef, nef,1);     # cortantes Qx y Qy

QQ = [ 
             16          20/3 (4*3^(1/2))/3            -16           20/3 -(4*3^(1/2))/3               0 8/3 0
 (16*3^(1/2))/3 (4*3^(1/2))/3             4 (16*3^(1/2))/3 -(4*3^(1/2))/3              4 -(32*3^(1/2))/3   0 8 ];
      
for e = 1:nef

    local x1, x2, x3, y1, y2, y3, Ae
    ## Calculo de la matriz de rigidez Ke
    x1 = xnod[LaG[e,1],X];              y1 = xnod[LaG[e,1],Y];
    x2 = xnod[LaG[e,2],X];              y2 = xnod[LaG[e,2],Y];
    x3 = xnod[LaG[e,3],X];              y3 = xnod[LaG[e,3],Y];

    idx = [ gdl[LaG[e,1],:]; gdl[LaG[e,2],:]; gdl[LaG[e,3],:] ]

    for i = 1:nGP
        # calculo de los momentos Mx, My y Mxy
        x = (1-L2[i]-L3[i])*x1 + L2[i]*x2 + L3[i]*x3;
        y = (1-L2[i]-L3[i])*y1 + L2[i]*y2 + L3[i]*y3;

        sigma_b[e,i] = -Dbe*L(x,y)*inv_A[e]*aa[idx]; 
    end
    # calculo de cortantes Qx y Qy: son constantes para todo el EF
    QxQy[e] = -D*QQ*aa[idx];
end

## Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1)  # número de elementos adyacentes
Mx  = zeros(nno,1)
My  = zeros(nno,1)
Mxy = zeros(nno,1)
Qx  = zeros(nno,1)
Qy  = zeros(nno,1)

# matriz de extrapolación

A = [ 9/4 5/4 -5/4 -5/4
       -9 5/2  5/2    5
       -9 5/2    5  5/2 ];

for e = 1:nef

    Mx[LaG[e,:],:] .+=     A * [sigma_b[e,1][1]
                                sigma_b[e,2][1]
                                sigma_b[e,3][1]
                                sigma_b[e,4][1] ]

    My[LaG[e,:],:] .+=     A * [sigma_b[e,1][2]
                                sigma_b[e,2][2]
                                sigma_b[e,3][2]
                                sigma_b[e,4][2] ]
                                            
    Mxy[LaG[e,:],:] .+=    A * [sigma_b[e,1][3]
                                sigma_b[e,2][3]
                                sigma_b[e,3][3]
                                sigma_b[e,4][3] ]

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
