#!/usr/bin/python3
from numpy import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def construir_matriz_newton(T):
    """Construye la matriz triangular inferior de Newton"""
    n = len(T)
    A = zeros((n, n))

    for i in range(n):
        A[i, 0] = 1  
        for j in range(1, n):
            if j > i:
                A[i, j] = 0
            else:
                prod = 1
                for k in range(j):
                    prod *= (T[i] - T[k]) 
                A[i, j] = prod
    return A

def sustitucion_adelante(A, b):
    """Resuelve Ax = b donde A es triangular inferior"""
    n = len(b)
    x = zeros(n)
    for i in range(n):
        suma = 0
        for j in range(i):
            suma += A[i, j] * x[j]
        x[i] = (b[i] - suma) / A[i, i]
    return x

T = array([500, 700, 900, 1100, 1300, 1500, 1700, 1900], dtype=float)
C = array([305, 301.6, 370.04, 379.08, 384.46, 395.5, 305.95, 307], dtype=float)

A = construir_matriz_newton(T)
b = C.copy()

print("\nMatriz de la base de polinomios de Newton:")
for fila in A:
    print("[", end=" ")
    for elemento in fila:
        print(f"{elemento:.2e}", end="  ")
    print("]")

coeficientes = sustitucion_adelante(A, b)

def evaluar_polinomio(t, T, coef):
    n = len(coef)
    resultado = coef[-1]
    for i in range(n-2, -1, -1):
        resultado = resultado * (t - T[i]) + coef[i] 
    return resultado

T_suave = linspace(min(T), max(T), 100)
C_interp = [evaluar_polinomio(t, T, coeficientes) for t in T_suave]

plt.figure(figsize=(10, 6))
plt.plot(T, C, 'ro', markersize=8, label='Datos originales')
plt.plot(T_suave, C_interp, 'b-', linewidth=2, label='Polinomio de Newton')
plt.title('Interpolación mediante Método de Newton - Lista 3', fontsize=14)
plt.ylabel('y', fontsize=12)
plt.xlabel('x', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

print("\nCoeficientes del polinomio interpolante:")
for i, c in enumerate(coeficientes):
    print(f"x_{i} = {c:.6e}")

plt.savefig('interpolacion_newton_matriz.pdf')
plt.show()









import numpy as np
import matplotlib.pyplot as plt

def construir_matriz_monomios(T):
    """Construye la matriz de monomios (Vandermonde)"""
    n = len(T)
    A = np.zeros((n, n))
    for i, t in enumerate(T):
        for j in range(n):
            A[i, j] = t**j
    return A

def evaluar_polinomio(t, coef):
    """Evalúa el polinomio en t dado un vector de coeficientes monomiales"""
    resultado = 0
    for j, a in enumerate(coef):
        resultado += a * (t**j)
    return resultado
# Datos
T = np.array([500, 700, 900, 1100, 1300, 1500, 1700, 1900], dtype=float)
C = np.array([305, 301.6, 370.04, 379.08, 384.46, 395.5, 305.95, 307], dtype=float)

# Construcción de A y resolución del sistema
A = construir_matriz_monomios(T)
coeficientes = np.linalg.solve(A, C)

# Evaluación del polinomio en una malla
T_suave = np.linspace(min(T), max(T), 200)
C_interp = [evaluar_polinomio(t, coeficientes) for t in T_suave]

# Gráfico
plt.figure(figsize=(10, 6))
plt.plot(T, C, 'ro', label='Datos originales')
plt.plot(T_suave, C_interp, 'b-', label='Polinomio base monomios')
plt.title('Interpolación con base monomial', fontsize=14)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.savefig('interpolacion_monomios.pdf')
plt.show()

# Mostrar coeficientes
print("\nCoeficientes del polinomio interpolante (base monomial):")
for i, c in enumerate(coeficientes):
    print(f"a_{i} = {c:.6e}")






    #!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

def lagrange_basis(x, i, X):
    """
    Calcula el polinomio base de Lagrange L_i(x).
    X: nodos (array con los x_j)
    i: índice del polinomio base
    """
    L = 1.0
    n = len(X)
    for j in range(n):
        if j != i:
            L *= (x - X[j]) / (X[i] - X[j])
    return L

def lagrange_interpol(x, X, Y):
    """
    Evalúa el polinomio de Lagrange en el punto x.
    X: nodos
    Y: valores de la función en los nodos
    """
    n = len(X)
    P = 0.0
    for i in range(n):
        P += Y[i] * lagrange_basis(x, i, X)
    return P

# Datos originales
T = np.array([500, 700, 900, 1100, 1300, 1500, 1700, 1900], dtype=float)
C = np.array([305, 301.6, 370.04, 379.08, 384.46, 395.5, 305.95, 307], dtype=float)

# Evaluación del polinomio en una malla
T_suave = np.linspace(min(T), max(T), 200)
C_interp = [lagrange_interpol(t, T, C) for t in T_suave]

# Gráfico
plt.figure(figsize=(10, 6))
plt.plot(T, C, 'ro', markersize=8, label='Datos originales')
plt.plot(T_suave, C_interp, 'b-', linewidth=2, label='Polinomio de Lagrange')
plt.title('Interpolación mediante Lagrange - Lista 3', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.savefig('interpolacion_lagrange.pdf')
plt.show()