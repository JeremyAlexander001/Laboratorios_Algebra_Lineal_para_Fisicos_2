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

