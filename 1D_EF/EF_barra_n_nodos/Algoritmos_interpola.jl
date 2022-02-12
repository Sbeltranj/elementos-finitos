using Interpolations, PyPlot, PyCall, MappedArrays

#https://juliamath.github.io/Interpolations.jl/latest/api/

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## Se crea un espacio para hacer clic y definir los nodos del EF
figure()
title("Haciendo clic con el ratón, defina los 10 puntos a interpolar")
axis([-5, 5, -5, 5])

coord = ginput(10)
xi = collect(1:1:10)
xnod = mappedarray(coord->first(coord),coord)
ynod = mappedarray(coord->last(coord),coord)


itp = interpolate(ynod, BSpline(Linear()))






