## Instala los paquetes de JULIA, utilizados en este repositorio 
## Santiago Beltrán Jaramillo

import Pkg

packages = ["Plots", "Polynomials", "PyPlot", "LinearAlgebra", "Statistics", "SparseArrays", "PyCall",
             "XLSX", "BoundaryValueDiffEq", "PlotlyJS", "MappedArrays", "WriteVTK", "Interpolations", "CubicSplines" ]

for e = 1:length(packages)
 Pkg.add(packages[e])
end