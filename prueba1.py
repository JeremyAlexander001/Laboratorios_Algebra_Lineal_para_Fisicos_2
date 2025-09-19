#! /usr/bin/python3
from numpy import *
import matplotlib.pyplot as plt # from pylab import plot,show
import warnings
warnings.filterwarnings("ignore")

#primeras 5 lineas, con condiciones del profesor

from numpy import zeros
datos = zeros([8,2], float)
from numpy import loadtxt
datos= loadtxt("data.txt", float)
print(datos)

#iniciamos con los datos, extraidos del archivo data.txt

x = datos[:,0]   # primera columna -> Q(l/h)
y = datos[:,1]   # segunda columna -> N(w)

#queremos hacer esta cosa p(t) = a0 + a1 t + a2 t^2 + ... + a7 t^7

n=8
#numero de datos (para despues generalizar)

A = [] #aqui se guardan los valores de la matriz A


for t in x: #t toma los valores de x, es decir, los valores de Q
    fila = []  #aqui se guardan los valores de cada fila
    for j in range(n): #j va de 0 a n-1
        fila.append(t**j)   # t^j
    A.append(fila) #agregamos la fila a la matriz A


#imprimimos la matriz A
for fila in A:  
    print(fila) 