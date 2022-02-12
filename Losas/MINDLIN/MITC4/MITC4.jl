## 
# Cálculo de los desplazamientos en una placa utilizando la teoría de
# Reissner-Mindlin y el elemento finito de placa MITC4
#
# Algoritmo documentado en:
# Katili, I., Batoz, J.-L., Maknun, J. and Lardeur, P. (2018), A comparative 
# formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements with new 
# numerical results based on s-norm tests. Computers & Structures, 204:
# 48-64. https://doi.org/10.1016/j.compstruc.2018.04.001
#
# y
#
# Bathe, K.-J. and Dvorkin, E.N. (1985), A four-node plate bending element 
# based on Mindlin/Reissner plate theory and a mixed interpolation. Int. J.
# Numer. Meth. Engng., 21: 367-383. https://doi.org/10.1002/nme.1620210213

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/losas/Mindlin/MITC4/EF_MITC4_Katili_et_al.m

# Programa elaborado en JULIA 1.7.1
# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

## cargamos paquetes:
using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK

include("losa.jl")  #para los gráficos
include("dib_DKMQ.jl")  #para los gráficos
close("all")          #cerrar ventanas

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

global nef = size(LaG,1)
## defino las variables/constantes
X = 1; Y = 2; Z = 3; # un par de constantes que ayudaran en la 
ww= 1; tx= 2; ty= 3; # lectura del código

E  = 210e9;       # [Pa]    módulo de elasticidad = 210GPa
nu = 0.3;         #         coeficiente de Poisson
h  = 0.05;        # [m]     espesor de la losa
q  = -10000;      # [N/m^2] carga
nno  = length(xnod[:,1]) #nno número de nodos

## definición de los grados de libertad
ngdl = 3*nno  
gdl  = [[1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]']    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,3)

## Se dibuja la malla de elementos finitos. 
figure(1)
cg = zeros(nef, 2) # almacena el centro de gravedad
nef   = size(LaG,1)
for e = 1:nef
    nod_ef = LaG[e, [1, 2, 3, 4, 1]]
    
    plt.plot(xnod[nod_ef,X], xnod[nod_ef,Y],
          color="k", linestyle="-")

    # Cálculo de la posición del centro de gravedad 
    cg[e,:] = [ mean(xnod[nod_ef,X]) mean(xnod[nod_ef,Y]) ]

    plt.text(cg[e, X], cg[e, Y], "$e", fontsize=5, color=[1,0,0],
        horizontalalignment="center", verticalalignment="center")

end

title("Malla de elementos finitos")
plt.plot(xnod[:,X], xnod[:,Y], "b.")


## Parámetros de la cuadratura de Gauss-Legendre
# se asumirá aquí el mismo orden de la cuadratura tanto en la dirección de
# xi como en la dirección de eta
include("gauss_legendre.jl")
n_gl = 2;                 # orden de la cuadratura de Gauss-Legendre
x_gl, w_gl = gausslegendre_quad(n_gl);

## Se leen las funciones de forma N y P y sus derivadas dN_dxi, dN_deta, 
#  dP_dxi, dP_deta
include("funciones_de_forma_Katili.jl");

## matrices constitutivas
Db = (E*h^3/(12*(1-nu^2)));   # plate rigidity
Hb = Db * [ 1  nu 0           # matriz constitutiva de flexión generalizada
            nu 1  0           # (Dbe en la nomenclatura del curso) 
            0  0  (1-nu)/2 ]; 

G  = E/(2*(1+nu));     # modulo de cortante
Hs = (5/6)*G*h*Matrix{Float64}(I, 2, 2); # matriz constitutiva de cortante generalizada (Dse)

## ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global
K   = spzeros(ngdl,ngdl)                 # matriz de rigidez global como RALA (sparse)
idx = Array{Array{Int64}}(undef, nef,1)
f   = zeros(ngdl,1);                     # vector de fuerzas nodales equivalentes global
N   = Array{Any}(undef,nef,n_gl, n_gl);

Bb = Array{Any}(undef,nef,n_gl,n_gl);  # matrices de deformación generalizada de flexión
Bs = Array{Any}(undef,nef,n_gl,n_gl); # matrices de deformación generalizada de cortante

for e = 1:nef               # ciclo sobre todos los elementos finitos
    ## Longitudes de los lados, cosenos y senos (Figura 4)
    xe = xnod[LaG[e,:],X];       ye = xnod[LaG[e,:],Y];
    x21 = xe[2] - xe[1];         y21 = ye[2] - ye[1]; 
    x32 = xe[3] - xe[2];         y32 = ye[3] - ye[2];
    x43 = xe[4] - xe[3];         y43 = ye[4] - ye[3];    
    x14 = xe[1] - xe[4];         y14 = ye[1] - ye[4];
    xji = [ x21 x32 x43 x14 ];   yji = [ y21 y32 y43 y14 ];   
    
    Lk = hypot.(xji, yji);      Ck =xji./Lk;      Sk = yji./Lk;
    
    ## Ciclo sobre los puntos de Gauss para calcular Kbe, Kse y fe
    Kbe = zeros(12,12);
    Kse = zeros(12, 12);
    fe  = zeros(12,1);
    det_Je = zeros(n_gl,n_gl); # almacenara los Jacobianos
    
    for pp = 1:n_gl
        for qq = 1:n_gl           
            ## Se evalúan las funciones de forma y sus derivadas en los 
            # puntos de Gauss
            xi_gl  = x_gl[pp];            eta_gl = x_gl[qq];

            NN       = Nforma(xi_gl, eta_gl);
            ddN_dxi  = dN_dxi(xi_gl, eta_gl);       
            ddN_deta = dN_deta(xi_gl, eta_gl);       

                                   
            ## Matriz jacobiana, su inversa y determinante
            # Se ensambla la matriz jacobiana
            dx_dxi  = sum(ddN_dxi .*xe);   dy_dxi  = sum(ddN_dxi .*ye);
            dx_deta = sum(ddN_deta.*xe);   dy_deta = sum(ddN_deta.*ye);
            
            Je = [ dx_dxi    dy_dxi
                   dx_deta   dy_deta ];
            
            # Se calcula su inversa
            inv_Je = inv(Je);
            j11 = inv_Je[1,1];              j12 = inv_Je[1,2];
            j21 = inv_Je[2,1];              j22 = inv_Je[2,2];                       
                
            # y su determinante (el Jacobiano)
            det_Je[pp,qq] = det(Je);

            # Se ensambla la matriz de funciones de forma N
            N[e,pp,qq] = zeros(3,12);
            for i = 1:4         
               N[e,pp,qq][:,[3*i-2 3*i-1 3*i]] = [ 
                   NN[i]    0           0
                   0       NN[i]        0
                   0       0           NN[i] ];
            end

            ## Se calcula la matriz de deformación por flexión Bb (ec. 40)
            Bb_ = zeros(3,12)
            for i = 1:4                
                dNi_dx = j11*ddN_dxi[i] + j12*ddN_deta[i]; # = ai
                dNi_dy = j21*ddN_dxi[i] + j22*ddN_deta[i]; # = bi               
                Bb_[:,[3*i-2 3*i-1 3*i]] = [ 0   dNi_dx        0
                                                     0        0   dNi_dy
                                                     0   dNi_dy   dNi_dx ];
            end

            Bb[e,pp,qq] = Bb_
            ## Se calcula la matriz de deformación por cortante Bs
            # Ecuación 52
            # Nota: esta ecuación se calculo en "demos_MITC4.m"

            Ng_Ag_Au = [        eta_gl/4 - 1/4       xi_gl/4 - 1/4
                         -(x21*(eta_gl - 1))/8 (x14*(xi_gl - 1))/8
                         -(y21*(eta_gl - 1))/8 (y14*(xi_gl - 1))/8
                                1/4 - eta_gl/4     -xi_gl/4 - 1/4
                         -(x21*(eta_gl - 1))/8 (x32*(xi_gl + 1))/8
                         -(y21*(eta_gl - 1))/8 (y32*(xi_gl + 1))/8
                                eta_gl/4 + 1/4       xi_gl/4 + 1/4
                         -(x43*(eta_gl + 1))/8 (x32*(xi_gl + 1))/8
                         -(y43*(eta_gl + 1))/8 (y32*(xi_gl + 1))/8
                              -eta_gl/4 - 1/4       1/4 - xi_gl/4
                         -(x43*(eta_gl + 1))/8 (x14*(xi_gl - 1))/8
                         -(y43*(eta_gl + 1))/8 (y14*(xi_gl - 1))/8 ]';
            
            Bs[e,pp,qq] = inv_Je*Ng_Ag_Au;
            
            ## se arma la matriz de rigidez del elemento e por flexión (eq. 45)
            Kbe +=  Bb[e,pp,qq]'*Hb*Bb[e,pp,qq]*det_Je[pp,qq]*w_gl[pp]*w_gl[qq];
            
            ## se arma la matriz de rigidez del elemento e por cortante (eq. 47)
            Kse += Bs[e,pp,qq]'*Hs*Bs[e,pp,qq]*det_Je[pp,qq]*w_gl[pp]*w_gl[qq];
            
            ## vector de fuerzas nodales equivalentes        
            if (xe[1] >= 0.9999 && xe[2] <= 1.501) && (ye[2] >= 0.9999 && ye[3] <= 2.001)
                fe += N[e,pp,qq]'*[q 0 0]'*det_Je[pp,qq]*w_gl[pp]*w_gl[qq];
            end
        end
    end
    
    ## ensamblaje matricial
    idx[e] = [ gdl[LaG[e,1],:]; gdl[LaG[e,2],:]; gdl[LaG[e,3],:]; gdl[LaG[e,4],:] ]
    K[idx[e],idx[e]] +=  Kbe + Kse;
    f[idx[e],:]      +=  fe;
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
a = zeros(ngdl);   a[c]  = ac;   a[d] = ad   # desplazamientos
q  = zeros(ngdl);  q[c]  = qd;   q[d] = qc   # fuerzas nodales equivalentes

## Dibujar deformada:
aa_ =  reshape(a,3,nno)'
a_  = aa_[:,X]
NL1, NL2, NL3, NL4 = 1,2,3,4

@pyimport matplotlib.tri as mtri
triangles = Vector{Vector{Int64}}(undef, 2*nef)

for e = 1:nef
    # se arma la matriz de correspondencia (LaG) de la malla
    triangles[2*e - 1] = LaG[e, [NL1, NL2, NL4]] .- 1
    triangles[2*e - 0] = LaG[e, [NL2, NL3, NL4]] .- 1
 end

triang = mtri.Triangulation(xnod[:,X], xnod[:,Y], triangles=triangles) 

fig = plt.figure()
esc = 0.8  #escala diagrama 
title("Estructura deformada $(esc) veces")
ax = fig.add_subplot(projection="3d")
ax.set_box_aspect((2, 4, esc)) 
ax.plot_trisurf(triang, a_, cmap="bwr")
plt.tight_layout() 

## Se calcula para cada elemento el vector de momentos en los puntos
## de Gauss (ecuacion 49)
MxMyMxy = Array{Any}(undef,nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            MxMyMxy[e,pp,qq] = Hb*Bb[e,pp,qq]*a[idx[e]];
        end
    end
end

## Se calcula para cada elemento el vector de cortantes en los puntos
## de Gauss 

QxQy = Array{Any}(undef,nef,n_gl,n_gl);
for e = 1:nef
    for pp = 1:n_gl
        for qq = 1:n_gl
            QxQy[e,pp,qq] = Hs*Bs[e,pp,qq]*a[idx[e]];
        end
    end
end

## Se extrapolan los momentos y cortantes a los nodos
num_elem_ady = zeros(nno,1);  # numero de elementos adyacentes
Mx  = zeros(nno,1);
My  = zeros(nno,1);
Mxy = zeros(nno,1);
Qx  = zeros(nno,1);
Qy  = zeros(nno,1);


A = [ 
   3^(1/2)/2 + 1            -1/2            -1/2   1 - 3^(1/2)/2
            -1/2   1 - 3^(1/2)/2   3^(1/2)/2 + 1            -1/2
   1 - 3^(1/2)/2            -1/2            -1/2   3^(1/2)/2 + 1
            -1/2   3^(1/2)/2 + 1   1 - 3^(1/2)/2            -1/2 ];

for e = 1:nef                             
    
    Mx[LaG[e,:],:] .+=   A *   [ MxMyMxy[e,1,1][1]
                                                MxMyMxy[e,1,2][1]
                                                MxMyMxy[e,2,1][1]
                                                MxMyMxy[e,2,2][1] ];

    My[LaG[e,:],:] .+=    A * [ MxMyMxy[e,1,1][2]
                                                MxMyMxy[e,1,2][2]
                                                MxMyMxy[e,2,1][2]
                                                MxMyMxy[e,2,2][2] ];
                                        
    Mxy[LaG[e,:],:] .+=   A * [ MxMyMxy[e,1,1][3]
                                                MxMyMxy[e,1,2][3]
                                                MxMyMxy[e,2,1][3]
                                                MxMyMxy[e,2,2][3] ];
                                            
    num_elem_ady[LaG[e,:],:] .+=  1;
end 


for e = 1:nef                             

    Qx[LaG[e,:],:] .+= A * [ QxQy[e,1,1][1]
                                                QxQy[e,1,2][1]
                                                QxQy[e,2,1][1]
                                                QxQy[e,2,2][1] ];

    Qy[LaG[e,:],:] .+= A * [ QxQy[e,1,1][2]
                                                QxQy[e,1,2][2]
                                                QxQy[e,2,1][2]
                                                QxQy[e,2,2][2] ];
end

## Alisado (promedio de los momentos y cortantes en los nodos)
Mx  =  Mx./num_elem_ady;  
My  =  My./num_elem_ady;  
Mxy = Mxy./num_elem_ady;   
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


## se dibujan los gráficos:
#Momentos Mx, My, Mxy  
figure(4)
subplot(131);plot_mom_Q_ang(xnod,[My], [],[L"Momento Mx(kN-m/m)"])
subplot(132); plot_mom_Q_ang(xnod,[Mx], [],[L"Momento My(kN-m/m)"])
subplot(133);plot_mom_Q_ang(xnod,[Mxy], [],[L"Momento Mxy(kN-m/m)"])

#Momentos principales
figure(5)
subplot(131);plot_mom_Q_ang(xnod,[Mf1_xy], [ang_],[L"Mf1_{xy}(kN-m/m)"])
subplot(132);plot_mom_Q_ang(xnod,[Mf2_xy], [ang_.+pi/2],[L"Mf2_{xy}(kN-m/m)"])
subplot(133);plot_mom_Q_ang(xnod,[Mt_max], [ang_.+pi/4, ang_.-pi/4],[L"Mt_{max}(kN-m/m)"])

#Cortantes Qx, Qy, Qmax 

figure(6)
subplot(131);plot_mom_Q_ang(xnod,[Qx], [],[L"Q_x(kN/m)"])
subplot(132);plot_mom_Q_ang(xnod,[Qy], [],[L"Q_y(kN/m)"])

## Se calculan y grafican los cortantes máximos, junto con su angulo de inclinacion

Q_max = hypot.(Qx, Qy);
ang   = atan.(Qy, Qx);
subplot(133);plot_mom_Q_ang(xnod,[Q_max], [ang],[ L"Q_{max}(kN/m)"])


#calculos wood_armer
include("wood_armer.jl")
wood = hcat(collect.(WoodArmer.(Mx, My, Mxy))...)'

#Diseño de wood y armer:
dibujar_wood_armer(xnod,[wood[:,1], wood[:,2], wood[:,3], wood[:,4]],
                [L"Momentos M_x^* sup", L"Momentos M_y^* sup",  L"Momentos M_x^* inf", L"Momentos M_y^* inf"]) 

## comparación solución analítica
u = 0.5; v = 1; xi = 1.25; eta = 1.5;
qdist = -10000;
err = zeros(nno,1);
MEF = zeros(nno,1);
analitica = zeros(nno,1);
ww = 1;

include("cal_w.jl")
for i = 1:nno
    MEF[i] = aa_[i,ww];
    analitica[i] = calc_w(xnod[i,X], xnod[i,Y], E, nu, h, 2, 4, qdist, u, v, xi, eta);
    err[i] = abs((MEF[i]-analitica[i])/analitica[i]);
end

println("Observe que al comparar ambos métodos los errores relativos máximos son:")
println(maximum(filter(!isnan,err)))
println("Es decir son extremadamente pequeños !!")
println()






