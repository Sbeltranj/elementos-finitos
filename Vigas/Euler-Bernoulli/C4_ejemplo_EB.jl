#JULIA 1.6.3
import XLSX
using SparseArrays
using PyPlot

# Programa para el cálculo de vigas de Euler-Bernoulli.

# %%Defino las constantes y variables
Y   = EF             = nodo       = 1
NL1 = TH             = direccion  = nodos = 2
NL2 = desplazamiento = fuerza_pun = 3
E_  = 4; I_ = 5; G_ = 6; Aas      = 7


#Nombre archivo EXCEL
filename = "viga_Uribe_Escamilla_ej_5_5.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(filename, "xnod")

# %%Se lee la posición de los nodos
T    = hcat(columns...)  
xnod = T[:,nodos]                       # Posición de los nodos
L    = diff(xnod)                      # Longitud de cada EF


nno  = length(xnod);                   # número de nodos
nef  = nno - 1;                        # número de elementos finitos (EF)
ngdl = 2*nno;                          # número de grados de libertad
gdl  = [ [1:2:ngdl]' [2:2:ngdl]' ];    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,2)


# %%Se leen la matriz de conectividad (LaG), el modulo de elasticidad, las 
# propiedades del material y las cargas
columns, labels = XLSX.readtable(filename, "LaG_EI_q")
T = hcat(columns...)

idxEF = T[:,EF]
LaG   = T[:,NL1:NL2]        # Definición de EFs respecto a nodos


E     = T[ :,E_];            # módulo de elasticidad E del EF
I     = T[ :,I_];            # momento de inercia Iz del EF
G     = T[ :,G_];            # módulo de rigidez (para viga de Timoshenko)
Aast  = T[ :,Aas];           # área de cortante (para viga de Timoshenko)
fz    = T[:,8:9];        # relación de las cargas distribuidas
fz = coalesce.(fz, 0.0)      # reemplazo los missing con ceros

#%%Relación de los apoyos

columns, labels = XLSX.readtable(filename, "restric")
T = hcat(columns...)

idxNODO = T[:,nodo]
dirdesp = T[:,direccion];
ac      = T[:,desplazamiento]; # desplazamientos conocidos
ac      = ac.*0.0

# %%Grados de libertad del desplazamiento conocidos y desconocidos
n_apoyos = length(idxNODO);  
c = zeros(n_apoyos, 1);        # GDL conocidos    

for i = 1:n_apoyos
  c[i,:] .= gdl[idxNODO[i], dirdesp[i]]
end

c = round.(Int, c)             # se convierte en Int64
c = vec(c)                     # ahora de matrix a vector
d =  setdiff(1:ngdl,c);        # GDL desconocidos

# %%Relación de cargas puntuales
columns, labels = XLSX.readtable(filename, "carga_punt")
T       = hcat(columns...)

idxNODO = T[:,nodo]
dirfp   = T[:,direccion];
fp      = T[:,fuerza_pun];; # desplazamientos conocidos


# %%Se colocan las fuerzas/momentos nodales en el vector de fuerzas nodales
# equivalentes global "f"
f_ini = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f_ini[gdl[idxNODO[i], dirfp[i]]] = fp[i];
end

# %%VIGA DE EULER-BERNOULLI:
# Con el programa "func_forma_euler_bernoulli.m" se calcularon:
#   Ke     = la matriz de rigidez de flexion del elemento e
#   fe     = el vector de fuerzas nodales equivalentes
#   Bb     = la matriz de deformaciones de flexion
#   N      = matriz de funciones de forma
#   dN_dxi = derivada de la matriz de funciones de forma con respecto a xi

# %%ensamblo la matriz de rigidez global y el vector de fuerzas nodales
# equivalentes global para la viga de Euler-Bernoulli


# %%Grados de libertad del desplazamiento conocidos y desconocidos
K_ini = spzeros(ngdl,ngdl); 
K = K_ini
f = f_ini
idx = Array{Array{Int64}}(undef, nef,1)      # grados de libertad del elemento e

for e = 1:nef  # Ciclo sobre todos los elementos finitos

    idx[e] = [gdl[LaG[e, 1], Y ]
              gdl[LaG[e, 1], TH]
              gdl[LaG[e, 2], Y ]
              gdl[LaG[e, 2], TH]]

    local    Le = L[e]

    # Matriz de rigidez de flexión del elemento e
    Ke = (E[e]*I[e]/Le^3) *[  12    6*Le   -12    6*Le    
                              6*Le  4*Le^2 -6*Le  2*Le^2 
                             -12   -6*Le    12   -6*Le    
                              6*Le  2*Le^2 -6*Le  4*Le^2 ]

    # vector de fuerzas nodales equivalentes de una carga trapezoidal 
    fe = [  (  Le*(7*fz[e,1] + 3*fz[e,2]))/20     # = Y1
            (Le^2*(3*fz[e,1] + 2*fz[e,2]))/60     # = M1
            (  Le*(3*fz[e,1] + 7*fz[e,2]))/20     # = Y2
           -(Le^2*(2*fz[e,1] + 3*fz[e,2]))/60 ]   # = M2
    # Se ensambla la matriz de rigidez K y el vector de fuerzas nodales
    # equivalentes f

    K[idx[e],idx[e]] += Ke
    f[idx[e]]        += fe # sumo a f global

end

# %%Se resuelve el sistema de ecuaciones
#%% extraigo las submatrices y especificó las cantidades conocidas
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K[c,c]; Kdc = K[d,c]; fc = f[d]
Kdd = K[d,d]; Kcd = K[c,d]; fd = f[c];

#%% resuelvo el sistema de ecuaciones

ad = Kdd\(fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad -fd

# armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);  # separó la memoria
a[c] = ac;          q[c] = qd;          
a[d] = ad;          #q[d] = qc = 0      


# %%Cálculo de los momentos flectores y las fuerzas cortantes
# M = se calcula en las raices del polinomio de GL de orden 2
# V = se calcula en el centro del EF (raiz del polinomio de GL de orden 1)
# se reserva la memoria
xmom = zeros(2,nef); # posición donde se calcula
mom  = zeros(2,nef); # momento flector
cor  = zeros(1,nef); # fuerza cortante

xi = [ -sqrt(1/3); sqrt(1/3) ]; # Raices del polinom de Legendre de grado 2

for e = 1:nef

   # longitud del elemento finito e
   local   Le = L[e];
   
   # matriz de deformaciones de flexión
   Bbe = [ (6 .*xi)/Le^2, (3 .*xi .- 1)/Le, -(6 .*xi)/Le^2, (3 .*xi .+ 1)/Le ]
   Bbe = hcat(Bbe...)
   
   # lugar donde se calcula el momento (centro del EF)
   xmom[:,e] = Le*xi'/2 .+ (xnod[LaG[e,1]] .+ xnod[LaG[e,2]])/2;
     
   # vector de desplazamientos nodales del elemento a^{(e)}
   local ae = a[idx[e]];
   
   mom[:,e] = E[e]*I[e] .*Bbe*ae;                 # momento flector   
   dN3_dxi3 = [ 3/2, (3*Le)/4, -3/2, (3*Le)/4 ]';
   cor[e]   = E[e]*I[e] *dN3_dxi3*(8/(Le^3))*ae; # fuerza cortante   

end


#%% se calculan los desplazamientos al interior de cada EF
nint = 10;           # número de puntos donde se interpolara dentro del EF
xi = collect(LinRange(-1,1,nint)); # coordenadas naturales


xx = Array{Array{Float64}}(undef, nef,1) # Interpol de posiciones (geometria) en el elemento
ww = Array{Array{Float64}}(undef, nef,1) # Interpol desplazamientos en el elemento
tt = Array{Array{Float64}}(undef, nef,1) # Interpol ángulo en el elemento

for e = 1:nef        # ciclo sobre todas los elementos finitos

   # longitud del elemento finito e
   local   Le = L[e]
   
   # Matriz de funciones de forma y su derivada
   local      N = [xi.^3/4 - (3*xi)/4 .+ 1/2,                   
                   -(Le*(- xi.^3/4 + xi.^2/4 + xi/4 .- 1/4))/2, 
                   - xi.^3/4 + (3*xi)/4 .+ 1/2,                 
                   -(Le*(- xi.^3/4 - xi.^2/4 + xi/4 .+ 1/4))/2] 
   N = hcat(N...)
    

   local dN_dxi = [(3*xi.^2)/4 .- 3/4,                          
                    -(Le*(- (3*xi.^2)/4 + xi/2 .+ 1/4))/2,      
                    3/4 .- (3*xi.^2)/4,                           
                    (Le*((3*xi.^2)/4 + xi/2 .- 1/4))/2 ]

   dN_dxi = hcat(dN_dxi...)

   #% vector de desplazamientos nodales del elemento a^{(e)}
   local   ae = a[idx[e]];

   # interpola sobre la geometria (coord naturales a geometricas)
   xx[e] = Le .*xi/2 .+ (xnod[LaG[e,1]] + xnod[LaG[e,2]])/2;
   
   # se calcula el desplazamiento al interior del elemento finito
   ww[e] = N*ae;
   
   # se calcula el ángulo al interior del elemento finito
   tt[e] = (atan.((dN_dxi*2/Le)*ae));
end



#%% Gráfico la solución analítica y la solución por el MEF
# %%1) grafico los desplazamientos de la viga

figure(1)                       # Lienzo nuevo
subplot(211)  

grid("on")                      # retícula
xlabel("Eje x [m]")             # Título del eje X
ylabel("Desplazamientos [m]")   # Título del eje Y
tight_layout()
title("Solución con el MEF para el desplazamiento ")


for e = 1:nef # ciclo sobre todos los elementos finitos
    defor = plot(xx[e], ww[e], color = :blue,
    label = nothing);           # gráfico solución por MEF
end

#%% 2) Gráfico de los ángulos de giro
subplot(212)
grid("on")
title("Solución con el MEF para el giro ", fontsize = 12)
ylabel("Giro (rad)")             # Título del eje Y
xlabel("Eje x [m]")              # Título del eje X

for e = 1:nef # ciclo sobre todos los elementos finitos
   giro =  plot(xx[e], tt[e], color = :blue,
   label = nothing);             # grafico solución por MEF
end

#%% 3) gráfico los momentos
figure(2)                        # cree un nuevo lienzo
subplot(212)

grid("on")                       # retícula
xlabel("Eje x [m]")              # Título del eje X
ylabel("Momento flector (kN-m)") # Título del eje Y
tight_layout()
title("Solución con el MEF para el momento flector")
plot(xmom[:], mom[:], color = :blue )


#%% 4) Gráfico de la fuerza cortante
subplot(211)

grid("on")
xlabel("Eje x [m]")             # Título del eje X
ylabel("Fuerza cortante (kN)")  # Título del eje Y
tight_layout()
title("Solución con el MEF para la fuerza cortante")

for e = 1:nef # ciclo sobre todos los elementos finitos
   plot([xnod[LaG[e,1]]; xnod[LaG[e,2]]],[cor[e]; cor[e]],
         color = :blue) #gráfico solución por MEF

end


display(figure(1))
display(figure(2))

gcf()

#%%Fin #falta agregar comparación por MAXIMA, e impresión bonita.
