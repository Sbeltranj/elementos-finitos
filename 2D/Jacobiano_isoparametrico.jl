# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

#adaptado y/o traducido:
#https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/2D/isoparametrico

# Tener instalada la librería de matplotlib (PYTHON-pip install matplotlib)
# además Pkg.add("PyPlot") en consola de JULIA
# https://discourse.julialang.org/t/recording-mouse-position-coordinates-in-plots/28848 

using PyPlot
using PyCall
using MappedArrays
close("all")          #cerrar ventanas
ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## Se crea un espacio para hacer clic y definir los nodos del EF
figure()
title("Haciendo clic con el ratón, defina los 16 nodos del EF")
axis([-5, 5, -5, 5])

coord = ginput(16)

xnod = mappedarray(coord->first(coord),coord)
ynod = mappedarray(coord->last(coord),coord)

## Espacio normalizado
fig = plt.figure()
title("Jacobiano Isoparamétrico")
gs = fig.add_gridspec(ncols=3, nrows=1)

#función meshgrid:
#tomado de: https://stackoverflow.com/questions/44581049/utilizing-ndgrid-meshgrid-functionality-in-julia

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

# se "discretiza" el espacio cada 0.05
delta = 0.05
xxi  = collect(-1:delta:1)
eeta = collect(-1:delta:1)
n    = length(xxi);
xi, eta = meshgrid(xxi, eeta)

# se grafica el EF en el espacio normalizado
ax1 = fig.add_subplot(gs[1,1])
for i = 1:n
    h1 = ax1.plot(xi[:, i], eta[:, i], "b")
    h2 = ax1.plot(xi[i, :], eta[i, :], "b")

    if i == 1 || i == (n)
        setp(h1, linewidth=4)
        setp(h2, linewidth=4)
    end
end
# se grafican los nodos en el espacio normalizado
xinod =  [-1, -1/3, 1/3,  1,    1,   1, 1, 1/3, -1/3, -1,  -1,   -1, -1/3,  1/3, 1/3, -1/3]
etanod = [-1,   -1,  -1, -1, -1/3, 1/3, 1,   1,    1,  1, 1/3, -1/3, -1/3, -1/3, 1/3,  1/3]
ax1.plot(xinod, etanod, "ro", markersize=12, linewidth=4)
ax1.axis([-1.1, 1.1, -1.1, 1.1])
ax1.set_title("Espacio normalizado")
ax1.set_aspect("equal", "box")

## Funciones de forma del elemento lagrangiano plano de 16 nodos (cuadrático)
#
# Numeración local:
#        ^ eta
#        |
#        |
# 10---9---8---7
#  |   |   |   |
# 11--16--15---6
#  |   |   |   |----> xi
# 12--13--14---5
#  |   |   |   |
#  1---2---3---4

Ni1 = (1/16).*(xi.-1).*(1 .-9 .*xi.^2);
Ni2 = (9/16).*(1 .-xi.^2).*(1 .-3 .*xi);
Ni3 = (9/16).*(1 .-xi.^2).*(1 .+3 .*xi);
Ni4 = (1/16).*(xi.+1).*(9 .*xi.^2 .-1);

Nj1 = (1/16).*(eta .-1).*(1 .-9 .*eta.^2);
Nj2 = (9/16).*(1 .-eta.^2).*(1 .-3 .*eta);
Nj3 = (9/16).*(1 .-eta.^2).*(1 .+3 .*eta);
Nj4 = (1/16).*(eta.+1).*(9 .*eta.^2 .-1);

N = Array{Array{Float64}}(undef,16)

N[10] = Ni1.*Nj4;   N[9]  = Ni2.*Nj4;   N[8]  = Ni3.*Nj4;   N[7]  = Ni4.*Nj4;
N[11] = Ni1.*Nj3;   N[16] = Ni2.*Nj3;   N[15] = Ni3.*Nj3;   N[6]  = Ni4.*Nj3;
N[12] = Ni1.*Nj2;   N[13] = Ni2.*Nj2;   N[14] = Ni3.*Nj2;   N[5]  = Ni4.*Nj2;
N[1]  = Ni1.*Nj1;   N[2]  = Ni2.*Nj1;   N[3]  = Ni3.*Nj1;   N[4]  = Ni4.*Nj1;

x = zeros(n,n)
y = zeros(n,n)
for i = 1:16
   x .+= N[i].*xnod[i];
   y .+= N[i].*ynod[i];
end

ax3 = fig.add_subplot(gs[1,2])
for i = 1:n
    h1 = ax3.plot(x[:, i], y[:, i], "b")
    h2 = ax3.plot(x[i, :], y[i, :], "b")
    if i == 0 || i == (n-1)
        plt.setp(h1, linewidth=4)
        plt.setp(h2, linewidth=4)
    end
end    
ax3.plot(xnod, ynod, "ro", markersize=12, linewidth=4)
ax3.set_aspect("equal", "box")


## derivadas de las funciones de forma
dN_dxi = zeros(n, n, 16)
dN_dxi[:,:,1]  = ((- 27 .*xi.^2 .+ 18 .*xi .+ 1).*(-9 .*eta.^3 .+ 9 .*eta.^2 + eta .- 1))/256;
dN_dxi[:,:,2]  = -(9 .*(-9 .*xi.^2 + 2 .*xi .+ 3).*(-9 .*eta.^3 .+ 9 .*eta.^2 + eta .- 1))/256;
dN_dxi[:,:,3]  = -(9 .*(9 .*xi.^2 .+ 2 .*xi .- 3).*(-9 .*eta.^3 .+ 9 .*eta.^2 + eta .- 1))/256;
dN_dxi[:,:,4]  = ((27 .*xi.^2 .+ 18 .*xi .- 1).*(-9 .*eta.^3 .+ 9 .*eta.^2 + eta .- 1))/256;
dN_dxi[:,:,5]  = (9 .*(3 .*eta .- 1).*(eta .- 1).*(eta .+ 1).*(27 .*xi.^2 .+ 18 .*xi .- 1))/256;
dN_dxi[:,:,6]  = (9 .*(27 .*xi.^2 .+ 18 .*xi .- 1).*(-3 .*eta.^3 - eta.^2 .+ 3 .*eta .+ 1))/256;
dN_dxi[:,:,7]  = -((27 .*xi.^2 .+ 18 .*xi .- 1).*(-9 .*eta.^3 .- 9 .*eta.^2 + eta .+ 1))/256;
dN_dxi[:,:,8]  = (9 .*(9 .*xi.^2 .+ 2 .*xi .- 3).*(-9 .*eta.^3 .- 9 .*eta.^2 + eta .+ 1))/256;
dN_dxi[:,:,9]  = (9 .*(- 9 .*xi.^2 .+ 2 .*xi .+ 3).*(-9 .*eta.^3 .- 9 .*eta.^2 + eta .+ 1))/256;
dN_dxi[:,:,10]  = -((-27 .*xi.^2 .+ 18 .*xi .+ 1).*(-9 .*eta.^3 .- 9 .*eta.^2 + eta .+ 1))/256;
dN_dxi[:,:,11]  = (9 .*(eta .- 1).*(3 .*eta .+ 1).*(eta .+ 1).*(27 .*xi.^2 .- 18 .*xi .- 1))/256;
dN_dxi[:,:,12]  = -(9 .*(eta .- 1).*(3 .*eta .- 1).*(eta .+ 1).*(27 .*xi.^2 .- 18 .*xi .- 1))/256;
dN_dxi[:,:,13]  = (81 .*(eta .- 1).*(3 .*eta .- 1).*(eta .+ 1).*(9 .*xi.^2 .- 2 .*xi .- 3))/256;
dN_dxi[:,:,14]  = (81 .*(9 .*xi.^2 .+ 2 .*xi .- 3).*(-3 .*eta.^3 + eta.^2 .+ 3 .*eta .- 1))/256;
dN_dxi[:,:,15]  = (81 .*(3 .*eta .+ 1).*(eta .- 1).*(eta .+ 1).*(9 .*xi.^2 .+ 2 .*xi .- 3))/256;
dN_dxi[:,:,16]  = -(81 .*(eta .- 1).*(3 .*eta .+ 1).*(eta .+ 1).*(9 .*xi.^2 - 2 .*xi .- 3))/256;


dN_deta = zeros(n, n, 16)
dN_deta[:,:,1] = ((-27 .*eta.^2 .+ 18 .*eta .+ 1).*(-9 .*xi.^3 .+ 9 .*xi.^2 + xi .- 1))/256;
dN_deta[:,:,2] = -(9 .*(xi .- 1).*(3 .*xi .- 1).*(xi .+ 1).*(27 .*eta.^2 .- 18 .*eta .- 1))/256;
dN_deta[:,:,3] = (9 .*(xi .- 1).*(3 .*xi .+ 1).*(xi .+ 1).*(27 .*eta.^2 .- 18 .*eta .- 1))/256;
dN_deta[:,:,4] = -((-27 .*eta.^2 .+ 18 .*eta .+ 1).*(-9 .*xi.^3 .- 9 .*xi.^2 + xi .+ 1))/256;
dN_deta[:,:,5] = (9 .*(-9 .*eta.^2 .+ 2 .*eta .+ 3).*(-9 .*xi.^3 .- 9 .*xi.^2 + xi .+ 1))/256;
dN_deta[:,:,6] = (9 .*(9 .*eta.^2 .+ 2 .*eta .- 3).*(-9 .*xi.^3 .- 9 .*xi.^2 + xi .+ 1))/256;
dN_deta[:,:,7] = -((27 .*eta.^2 .+ 18 .*eta .- 1).*(-9 .*xi.^3 .- 9 .*xi.^2 + xi .+ 1))/256;
dN_deta[:,:,8] = (9 .*(27 .*eta.^2 .+ 18 .*eta .- 1).*(-3 .*xi.^3 - xi.^2 .+ 3 .*xi .+ 1))/256;
dN_deta[:,:,9] = (9 .*(3 .*xi .- 1).*(xi .- 1).*(xi .+ 1).*(27 .*eta.^2 .+ 18 .*eta .- 1))/256;
dN_deta[:,:,10] = ((27 .*eta.^2 .+ 18 .*eta .- 1).*(-9 .*xi.^3 .+ 9 .*xi.^2 + xi .- 1))/256;
dN_deta[:,:,11] = -(9 .*(9 .*eta.^2 .+ 2 .*eta .- 3).*(-9 .*xi.^3 .+ 9 .*xi.^2 + xi .- 1))/256;
dN_deta[:,:,12] = -(9 .*(-9 .*eta.^2 .+ 2 .*eta .+ 3).*(-9 .*xi.^3 .+ 9 .*xi.^2 + xi .- 1))/256;
dN_deta[:,:,13] = (81 .*(xi .- 1).*(3 .*xi .- 1).*(xi .+ 1).*(9 .*eta.^2 .- 2 .*eta .- 3))/256;
dN_deta[:,:,14] = -(81 .*(xi .- 1).*(3 .*xi .+ 1).*(xi .+ 1).*(9 .*eta.^2 .- 2 .*eta .- 3))/256;
dN_deta[:,:,15] = (81 .*(3 .*xi .+ 1).*(xi .- 1).*(xi .+ 1).*(9 .*eta.^2 .+ 2 .*eta .- 3))/256;
dN_deta[:,:,16] = (81 .*(9 .*eta.^2 .+ 2 .*eta .- 3).*(-3 .*xi.^3 + xi.^2 .+ 3 .*xi .- 1))/256;

## Estas derivadas se calcularon con el siguiente código de MATLAB:

#=
Ni1 = (1/16)*(xi - 1)*(1 - 9*xi**2)
Ni2 = (9/16)*(1 - xi**2)*(1 - 3*xi)
Ni3 = (9/16)*(1 - xi**2)*(1 + 3*xi)
Ni4 = (1/16)*(xi + 1)*(9*xi**2 - 1)
Nj1 = (1/16)*(eta - 1)*(1 - 9*eta**2)
Nj2 = (9/16)*(1 - eta ** 2)*(1 - 3*eta)
Nj3 = (9/16)*(1 - eta ** 2)*(1 + 3*eta)
Nj4 = (1/16)*(eta + 1)*(9*eta**2 - 1)
N = 16 * [None]
N[ 9] = Ni1 * Nj4; N[ 8] = Ni2 * Nj4; N[ 7] = Ni3 * Nj4; N[6] = Ni4 * Nj4; 
N[10] = Ni1 * Nj3; N[15] = Ni2 * Nj3; N[14] = Ni3 * Nj3; N[5] = Ni4 * Nj3;
N[11] = Ni1 * Nj2; N[12] = Ni2 * Nj2; N[13] = Ni3 * Nj2; N[4] = Ni4 * Nj2;
N[ 0] = Ni1 * Nj1; N[ 1] = Ni2 * Nj1; N[ 2] = Ni3 * Nj1; N[3] = Ni4 * Nj1;
dN_dxi  = 16 * [None]
dN_deta = 16 * [None]
for i in range(16):
    dN_dxi[i] = simplify(sp.diff(N[i], xi))
    print(f'dN_dxi{i} = {char(dN_dxi[i])}')
for i in range(16):
    dN_deta[i] = simplify(sp.diff(N[i], eta))
    print(f'dN_dxi{i} = {char(dN_deta[i])}')
=#

## Calculo el determinante del Jacobiano
dx_dxi_  = zeros(n, n)
dy_dxi_  = zeros(n, n)
dx_deta_ = zeros(n, n)
dy_deta_ = zeros(n, n)

for i = 1:16
   global dx_dxi_  += dN_dxi[:,:,i].*xnod[i];
   global dy_dxi_  += dN_dxi[:,:,i].*ynod[i];   
   global dx_deta_ += dN_deta[:,:,i].*xnod[i];
   global dy_deta_ += dN_deta[:,:,i].*ynod[i];   
end

# Calculo el determinante del Jacobiano
# J = [ dx_dxi   dy_dxi
#       dx_deta  dy_deta ]
detJ = dx_dxi_.*dy_deta_ - dx_deta_.*dy_dxi_

# Se calcula el Jacobian ratio
if maximum(detJ) > 0
    JR = maximum(detJ) / minimum(detJ)
else
    JR = minimum(detJ) / maximum(detJ)
end


## se grafica el jacobiano
ax2 = fig.add_subplot(gs[1,3])
h = ax2.pcolor(xi, eta, detJ, cmap="jet")
tickscb = LinRange(minimum(detJ), maximum(detJ), 10)
fig.colorbar(h, shrink=0.5), 
ax2.set_title("Determinante de J")

if JR < 0 || JR > 40
    ax3.set_title("Esta forma no es adecuada para un EF.")
else
    ax3.set_title("Jacobian ratio = $(round(JR, digits = 2))")  
end

ax2.set_aspect("equal", "box")
ax2.contour(xi, eta, detJ, levels= [0] , linewidths=4)


if JR > 0
    println("Jacobian ratio (JR) = ",JR)
    println("min(detJ) = ",minimum(detJ))
    println("max(detJ) = ",maximum(detJ))
end

if maximum(detJ) < 0
    println("El EF está numerado en sentido contrario al especificado")
end
if JR < 0 || JR > 40
    println("Esta forma no es adecuada para un EF.")
end
##fin