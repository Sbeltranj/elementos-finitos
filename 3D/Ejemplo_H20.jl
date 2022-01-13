using LinearAlgebra, Statistics, SparseArrays,  WriteVTK
import XLSX


#Algunas variables, para facilitar lectura del código
X = nodo = 1
x_ = NL1 = direccion = Ee = direccion = Y = 2
fpuntual = desplazamiento = Z =3
z_  = 4
NL20 = 21
g = 9.81 #m/s2 

#Nombre archivo EXCEL
#nombre_archivo = "malla_H20_conexion.xlsx"
nombre_archivo = "malla_H20_viga.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(nombre_archivo, "xnod")

##posición de los nodos:
##Se lee la posición de los nodos
T    = hcat(columns...)  

# xnod: fila=número del nodo, columna=coordenada X_=1 o Y_=2
xnod = T[:,x_:z_] 
xnod = xnod.*1.0  # convert Type real to Float64

nno  = length(xnod[:,1])

## definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
columns, labels = XLSX.readtable(nombre_archivo, "LaG_mat")
T = hcat(columns...)

idx_LaG = [2, 10, 3, 11, 4, 12, 5, 13, 14, 15, 16, 17, 6, 18, 7, 19, 8, 20, 9, 21]
LaG   = T[:,idx_LaG]        # Definición de EFs respecto a nodos
nef   = size(LaG,1)
nnpe  = size(LaG,2)          # número de nodos por EF (=20)

if nnpe != 20
    error("Este código SOLO es sirve para EFs hexaédricos de 20 nodos (H20)")
end

## definición de los grados de libertad
ngdl = 3*nno            # número de grados de libertad por nodo = [X, Y, Z]
gdl  = [[1:3:ngdl]' [2:3:ngdl]' [3:3:ngdl]']    # grados de libertad
gdl  = reshape(hcat(gdl...)',nno,3)

## material
columns, labels = XLSX.readtable(nombre_archivo, "prop_mat")
T       = hcat(columns...)

E          = T[:,2]   # módulo de elasticidad E
nu         = T[:,3]  # coeficiente de Poisson
rhoe       = T[:,4]  # densidad del material
nmat       = length(E)
## Relación de cargas puntuales

columns, labels = XLSX.readtable(nombre_archivo, "carga_punt")
T  = hcat(columns...)

ncp     = size(T,1)        # número de cargas puntuales
idxNODO = T[:,2]        
dirfp   = T[:,direccion];  # dirección carga
fp      = T[:,fpuntual];   # vector de fuerzas nodales equivalentes global

f = zeros(ngdl,1);   # vector de fuerzas nodales equivalentes global

for i = 1:length(idxNODO)
   f[gdl[idxNODO[i], dirfp[i]]] = fp[i]
end

include("funcion_H20.jl")
## Cuadratura de Gauss-Legendre
# NOTA: se asumirá aquí el mismo orden de la cuadratura tanto en la dirección
#       de xi como en la dirección de eta
x_gl, w_gl = gausslegendre_quad_hexa(2) # 2x2x2
n_gl = length(w_gl)

### Funciones de forma (serendípitas) y sus derivadas del elemento rectangular
#   de 8 nodos:
Nforma   = N_H20
dN_dxi   = dN_dxi_H20
dN_deta  = dN_deta_H20
dN_dzeta = dN_dzeta_H20

# se hace la iniciación de la matriz de rigidez global y los espacios en memoria que
#  almacenarán las matrices de forma y de deformación

K   = spzeros(ngdl,ngdl)  # matriz de rigidez global como RALA (sparse)
N = Array{Any}(undef,nef,n_gl)      # contenedor para las matrices de forma
B = Array{Any}(undef,nef,n_gl)      # contenedor para las matrices de deformación

# matriz constitutiva del elemento para TENSIÓN PLANA
De = Array{Any}(undef,nmat,1)
be = Array{Any}(undef,nmat,1)

for i  = 1:nmat
    
    local d1, d2
    d1 = (1-nu[i])/(1-2*nu[i])
    d2 = nu[i]/(1-2*nu[i])

    global De
    De = E[i]/(1+nu[i]).*     [d1 d2 d2  0   0   0  
                               d2 d1 d2  0   0   0  
                               d2 d2 d1  0   0   0  
                                0  0  0  1/2 0   0  
                                0  0  0  0   1/2 0  
                                0  0  0  0   0   1/2]
    global be
    be = [0; -rhoe[i]*g; 0]  # [kgf/m³] vector de fuerzas másicas

end

for e = 1:nef          # ciclo sobre todos los elementos finitos
    # Calculo la matriz de rigidez y el vector de fuerzas nodales
    # equivalentes del elemento

    local Ke, fe, det_Je
    Ke = zeros(3*nnpe, 3*nnpe );
    fe = zeros(3*nnpe, 1);
    det_Je = zeros(n_gl, 1); # en esta matriz se almacenaran los Jacobianos
    
    for p = 1:n_gl

       xi_gl   = x_gl[p,1];
       eta_gl  = x_gl[p,2];
       zeta_gl = x_gl[p,3];
       
       # Se evalúan las funciones de forma en los puntos de integración
       # de Gauss-Legendre
       NNforma = Nforma(xi_gl, eta_gl, zeta_gl);
       
       # Se evalúan las derivadas de las funciones de forma en los puntos
       # de integración de Gauss-Legendre
       local xi_gl, eta_gl, NNforma, ddN_dxi, ddN_deta, xe, ye, ze, zeta_gl, ddN_dzeta

       xe = xnod[LaG[e,:],X]; ye = xnod[LaG[e,:],Y]; ze = xnod[LaG[e,:],Z];

       ddN_dxi   = dN_dxi(xi_gl, eta_gl, zeta_gl);    
       ddN_deta  = dN_deta(xi_gl, eta_gl, zeta_gl);   
       ddN_dzeta = dN_dzeta(xi_gl, eta_gl, zeta_gl);  

       local dx_dxi, dx_deta, dy_dxi, dy_deta, Je, dz_dxi, dz_deta, dx_dzeta, dy_dzeta, dz_dzeta
       dx_dxi   = sum(ddN_dxi  .*xe);   dy_dxi   = sum(ddN_dxi  .*ye);   dz_dxi   = sum(ddN_dxi  .*ze);
       dx_deta  = sum(ddN_deta .*xe);   dy_deta  = sum(ddN_deta .*ye);   dz_deta  = sum(ddN_deta .*ze);
       dx_dzeta = sum(ddN_dzeta.*xe);   dy_dzeta = sum(ddN_dzeta.*ye);   dz_dzeta = sum(ddN_dzeta.*ze);
       
       # Se ensambla la matriz Jacobiana del elemento
       Je = [ dx_dxi    dy_dxi    dz_dxi
              dx_deta   dy_deta   dz_deta
              dx_dzeta  dy_dzeta  dz_dzeta ];
       
       # Se calcula el determinante del Jacobiano
       det_Je[p] = det(Je);
       
       
       # las matrices de forma y de deformación se evalúan y se ensamblan
       # en el punto de Gauss

       N[e,p] = zeros(3,3*nnpe) #Array{Any}(undef,3,3*nnpe)
       B[e,p] = zeros(6,3*nnpe)

       # se ensamblan la matriz de rigidez del elemento y el vector de
       # fuerzas nodales equivalentes del elemento
       
       for i = 1:nnpe
          # Se ensambla la matriz de funciones de forma N
          N[e,p][:,[3*i-2 3*i-1 3*i]] = [ 
             NNforma[i]  0           0
             0           NNforma[i]  0
             0           0           NNforma[i] ];

             # Se ensambla la matriz de deformaciÓn del elemento B
             tmp = Je\[ ddN_dxi[i];  ddN_deta[i];  ddN_dzeta[i] ];

             local dNi_dx, dNi_dy, dNi_dz

             dNi_dx = tmp[1];
             dNi_dy = tmp[2];
             dNi_dz = tmp[3];

          B[e,p][:,[3*i-2 3*i-1 3*i]] = [ 
                dNi_dx 0      0        # aquí se ensambla
                0      dNi_dy 0        # y asigna la matriz
                0      0      dNi_dz   # B_i
                dNi_dy dNi_dx 0
                dNi_dz 0      dNi_dx
                0      dNi_dz dNi_dy ];
       end  

        # se arma la matriz de rigidez del elemento e
        Ke += B[e,p]'*De*B[e,p]*det_Je[p]*w_gl[p]
        #Y el vector local de fuerzas equivalentes
        fe += N[e,p]'*be*det_Je[p]*w_gl[p]
    end
    # se determina si hay puntos con jacobiano negativo, en caso tal se termina
    # el programa y se reporta
#=     if det_Je .< 0
        error("Hay puntos con det_Je negativo en el elemento ")
    end =#
    #idx: se asignan los grados de libertad según LaG
    local idx
    idx = Array{Int64}(undef, 1, 3*nnpe)

    for i = 1:nnpe
        idx[3*i-2:3*i] = gdl[LaG[e,i],:]
    end
    idx = vec(idx)
    
    #Se realiza el ensamblaje en el vector global f y la matriz K.
    f[idx,:]   += fe;
    K[idx,idx] += Ke;
    
end
 
## Muestro la configuración de la matriz K (K es rala)
#= figure(1)
spy(K)
title("Los puntos representan los elementos diferentes de cero") =#


## se relacionan las restricciones
columns, labels = XLSX.readtable(nombre_archivo, "restric")
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

## Se calcula para cada elemento las deformaciones y los esfuerzos
def = Array{Any}(undef,nef,n_gl)
esf = Array{Any}(undef,nef,n_gl)

for e = 1:nef

   local idx
   idx = Array{Int64}(undef, 1, 3*nnpe)

   for i = 1:nnpe
       idx[3*i-2:3*i] = gdl[LaG[e,i],:]
   end

   idx = vec(idx)
   ae = a[idx]            # desplazamientos de los gdl del elemento e
   
   for p = 1:n_gl
      def[e,p] = B[e,p]*ae     # calculo las deformaciones
      esf[e,p] = De*def[e,p]   # calculo los esfuerzos
   end
   
end


#matriz_extrapolacion_esfuerzos_H20
A = matriz_extrapolacion_esfuerzos_H20()

## Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno)  # número de elementos adyacentes
sx  = zeros(nno);   ex  = zeros(nno)
sy  = zeros(nno);   ey  = zeros(nno)
sz  = zeros(nno);   gxy = zeros(nno)
txy = zeros(nno);   gxz = zeros(nno)
txz = zeros(nno);   gyz = zeros(nno)
tyz = zeros(nno);   ez  = zeros(nno)


# se hace la extrapolación de los esfuerzos y las deformaciones en cada elemento
# a partir de las lecturas en los puntos de Gauss
## Me falta optimizar esta parte, con un ciclo for...

for e = 1:nef
   sx[LaG[e,:],:] .+=    A *  [ esf[e,1][1], esf[e,2][1], esf[e,3][1], esf[e,4][1], esf[e,5][1], esf[e,6][1], 
                               esf[e,7][1], esf[e,8][1] ]
   sy[LaG[e,:],:] .+=    A *  [ esf[e,1][2], esf[e,2][2], esf[e,3][2], esf[e,4][2], esf[e,5][2], esf[e,6][2], 
                               esf[e,7][2], esf[e,8][2] ]

   sz[LaG[e,:],:] .+=    A *  [ esf[e,1][3], esf[e,2][3], esf[e,3][3], esf[e,4][3], esf[e,5][3], esf[e,6][3], 
                               esf[e,7][3], esf[e,8][3] ]

   txy[LaG[e,:],:] .+=   A *  [ esf[e,1][4], esf[e,2][4], esf[e,3][4], esf[e,4][4], esf[e,5][4], esf[e,6][4], 
                               esf[e,7][4], esf[e,8][4] ]

   txz[LaG[e,:],:] .+=   A *  [ esf[e,1][5], esf[e,2][5], esf[e,3][5], esf[e,4][5], esf[e,5][5], esf[e,6][5], 
                               esf[e,7][5], esf[e,8][5] ]

   tyz[LaG[e,:],:] .+=   A *  [ esf[e,1][6], esf[e,2][6], esf[e,3][6], esf[e,4][6], esf[e,5][6], esf[e,6][6], 
                               esf[e,7][6], esf[e,8][6] ]
                                
   ex[LaG[e,:],:] .+=    A *  [ def[e,1][1], def[e,2][1], def[e,3][1], def[e,4][1], def[e,5][1], def[e,6][1], 
                               def[e,7][1], def[e,8][1] ]

   ey[LaG[e,:],:] .+=    A *  [ def[e,1][2], def[e,2][2], def[e,3][2], def[e,4][2], def[e,5][2], def[e,6][2], 
                               def[e,7][2], def[e,8][2] ]

   ez[LaG[e,:],:] .+=    A *  [ def[e,1][3], def[e,2][3], def[e,3][3], def[e,4][3], def[e,5][3], def[e,6][3], 
                               def[e,7][3], def[e,8][3] ]

   gxy[LaG[e,:],:] .+=   A *  [ def[e,1][4], def[e,2][4], def[e,3][4], def[e,4][4], def[e,5][4], def[e,6][4], 
                               def[e,7][4], def[e,8][4] ]

   gxz[LaG[e,:],:] .+=   A *  [ def[e,1][5], def[e,2][5], def[e,3][5], def[e,4][5], def[e,5][5], def[e,6][5], 
                               def[e,7][5], def[e,8][5] ]

   gyz[LaG[e,:],:] .+=   A *  [ def[e,1][6], def[e,2][6], def[e,3][6], def[e,4][6], def[e,5][6], def[e,6][6], 
                               def[e,7][6], def[e,8][6] ]

   num_elem_ady[LaG[e,:],:] .+= 1
end

## Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;        ex  =  ex./num_elem_ady;
sy  =  sy./num_elem_ady;        ey  =  ey./num_elem_ady;
sz  =  sz./num_elem_ady;        ez  =  ez./num_elem_ady;
txy = txy./num_elem_ady;        gxy = gxy./num_elem_ady;
txz = txz./num_elem_ady;        gxz = gxz./num_elem_ady;
tyz = tyz./num_elem_ady;        gyz = gyz./num_elem_ady;

## Se calculan para cada nodo los esfuerzos principales y sus direcciones
# %% Se calculan para cada nodo los esfuerzos principales y sus direcciones
s1 = zeros(nno);  n1 = zeros((nno, 3))
s2 = zeros(nno);  n2 = zeros((nno, 3))
s3 = zeros(nno);  n3 = zeros((nno, 3))

for i = 1:nno

    local esfppales, dirppales
    esfppales, dirppales = eigen(
                              [sx[i]   txy[i]    txz[i]  # matriz de esfuerzos
                               txy[i]   sy[i]    tyz[i]     # de Cauchy para 
                               txz[i]  tyz[i]    sz[i] ]) # theta = grados

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

## Export results ParaView EF_H20
cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, vec(LaG[e,[1 3 5 7 13 15 17 19] ]) ) for e = 1:nef]

vtkfile = vtk_grid("H20_element", xnod[:,X],xnod[:,Y], xnod[:,Z], cells) 

vtkfile["s_x"] = sx;    vtkfile["txy"] = txy;      vtkfile["ex"] = ex; vtkfile["uvw"] = a  ; vtkfile["gxz"] = gxz
vtkfile["s_y"] = sy;    vtkfile["txz"] = txz;      vtkfile["ey"] = ey; vtkfile["sv"]  = sv ; vtkfile["gyz"] = gyz
vtkfile["s_z"]  =sz;    vtkfile["tyz"] = tyz;      vtkfile["ez"] = ez; vtkfile["gxy"] = gxy; vtkfile["Tmax"] = tmax
outfiles = vtk_save(vtkfile)
println("Se han reportado los resultados en formato.vtu para ser visualizados en ParaView")