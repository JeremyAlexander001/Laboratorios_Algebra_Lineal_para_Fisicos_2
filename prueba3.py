#!/usr/bin/python3

from numpy import *

import numpy as np

import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore")

# -DATOS
T = np.array([500, 700, 900, 1100, 1300, 1500, 1700, 1900], dtype=float)
C = np.array([365, 361.6, 370.64, 379.08, 384.46, 395.5, 395.95, 397], dtype=float)

# ---MONOMIOS---

# construye la matriz de vandermonde
def construir_matriz_monomios(T):
    n = len(T)
    A = np.zeros((n, n))
    for i, t in enumerate(T):
        for j in range(n):
            A[i, j] = t**j
    return A

# EVALUAR EL POLINOMIO EN t dado un arreglo de coeficientes
def evaluar_polinomio_monico(t, coef):
    resultado = 0
    for j, a in enumerate(coef):
        resultado += a * (t**j)
    return resultado

# ---NEWTON---

#  MATRIZ TRIANGULAR INFERIOR DE NEWTON
def construir_matriz_newton(T):
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

# Sustitucion directa (forward)
def sustitucion_adelante(A, b):
    n = len(b)
    x = zeros(n)
    for i in range(n):
        suma = 0
        for j in range(i):
            suma += A[i, j] * x[j]
        x[i] = (b[i] - suma) / A[i, i]
    return x

def evaluar_polinomio_newton(t, T, coef):
    n = len(coef)
    resultado = coef[-1]
    for i in range(n-2, -1, -1):
        resultado = resultado * (t - T[i]) + coef[i] 
    return resultado

# ---LAGRANGE--- 

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

# ---- Monomios ----
A_vand = construir_matriz_monomios(T)
# Mostrando la matriz de vadermonde
print("--Matriz de vadermonde--")

for fila in A_vand:
    print("[", end=" ")
    for elemento in fila:
        print(f"{elemento:.2e}", end="  ")
    print("]")
    
coef_monomial = np.linalg.solve(A_vand, C)

print("\nCoeficientes del polinomio interpolante:")
for i, c in enumerate(coef_monomial):
    print(f"x_{i} = {c:.6e}")

# ---NEWTON---

A_Newton = construir_matriz_newton(T)
b = C.copy()
coef_newton = sustitucion_adelante(A_Newton, C.copy())
# MOSTRANDO LA MATRIZ DE NEWTON
print("\nMatriz de la base de polinomios de Newton:")
for fila in A_Newton:
    print("[", end=" ")
    for elemento in fila:
        print(f"{elemento:.2e}", end="  ")
    print("]")

## ---LAGRANGE---




T_suave = np.linspace(min(T), max(T), 200)

C_newton_interp = [evaluar_polinomio_newton(t, T, coef_newton) for t in T_suave]

C_monomial_interp = [evaluar_polinomio_monico(t, coef_monomial) for t in T_suave]

C_lagrange_interp = [lagrange_interpol(t, T, C) for t in T_suave]

# -- Gráfica comparativa --
plt.figure(figsize=(10, 6))
plt.plot(T, C, 'ro', markersize=8, label='Datos originales')
plt.plot(T_suave, C_newton_interp,'+', linewidth=2, label='Polinomio de Newton')
plt.plot(T_suave, C_monomial_interp, '--', linewidth=1.5, label='Polinomio base monomios')
plt.plot(T_suave, C_lagrange_interp, ':', linewidth=1.5, label='Polinomio de Lagrange')
plt.title('Comparación: Newton, Monomial y Lagrange', fontsize=14)
plt.ylabel('N(w)', fontsize=12)
plt.xlabel('Q(l/h)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.savefig('interpolaciones_fusionadas.pdf')
plt.show()
