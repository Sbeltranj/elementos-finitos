using Polynomials,LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK
import XLSX

## exportando geometría a Paraview
X, Y, Z = 1,2,3
x_  = 2
z_  = 4

#Nombre archivo EXCEL
nombre_archivo = "malla_H20_conexion.xlsx"

#se carga el libro.xlsx, con el nombre de la hoja "xnod"
columns, labels = XLSX.readtable(nombre_archivo, "xnod")

##posición de los nodos:
##Se lee la posición de los nodos
T    = hcat(columns...)  

# xnod: fila=número del nodo, columna=coordenada X_=1 o Y_=2
xnod = T[:,x_:z_] 
xnod = xnod.*1.0

## definición de elementos finitos con respecto a nodos
# LaG: fila=número del elemento, columna=número del nodo local
columns, labels = XLSX.readtable(nombre_archivo, "LaG_mat")
T = hcat(columns...)

LaG   = T[:,2:21]        # Definición de EFs respecto a nodos
nef   = size(LaG,1)

cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, vec(LaG[e,[1 2 3 4 5 6 7 8] ]) ) for e = 1:nef]

vtkfile = vtk_grid("H20_element", xnod[:,X],xnod[:,Y], xnod[:,Z], cells) 
outfiles = vtk_save(vtkfile)