#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Construcción de la matriz de Newton
def construir_matriz_newton(T):
    """Construye la matriz triangular inferior de Newton"""
    n = len(T)
    A = np.zeros((n, n))

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
    x = np.zeros(n)
    for i in range(n):
        suma = 0
        for j in range(i):
            suma += A[i, j] * x[j]
        x[i] = (b[i] - suma) / A[i, i]
    return x

# Datos originales
T = np.array([500, 700, 900, 1100, 1300, 1500, 1700, 1900], dtype=float)
C = np.array([365, 301.6, 370.04, 379.08, 384.46, 395.5, 305.95, 307], dtype=float)

# Malla común para graficar
T_suave = np.linspace(min(T), max(T), 200)

# Newton
A = construir_matriz_newton(T)
b = C.copy()

print("\nMatriz de la base de polinomios de Newton:")
for fila in A:
    print("[", end=" ")
    for elemento in fila:
        print(f"{elemento:.2e}", end="  ")
    print("]")

coeficientes_newton = sustitucion_adelante(A, b)

def evaluar_polinomio_newton(t, T, coef):
    n = len(coef)
    resultado = coef[-1]
    for i in range(n-2, -1, -1):
        resultado = resultado * (t - T[i]) + coef[i] 
    return resultado

print("\nCoeficientes del polinomio interpolante (Newton):")
for i, c in enumerate(coeficientes_newton):
    print(f"x_{i} = {c:.6e}")

C_interp_newton = [evaluar_polinomio_newton(t, T, coeficientes_newton) for t in T_suave]

# Monomial (Vandermonde)
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

A = construir_matriz_monomios(T)
coeficientes_monomial = np.linalg.solve(A, C)

C_interp_monomial = [evaluar_polinomio(t, coeficientes_monomial) for t in T_suave]

np.set_printoptions(precision=3, suppress=True, linewidth=150)
print("\nMatriz de monomios (Vandermonde):\n")
print(A)

print("\nCoeficientes del polinomio interpolante (base monomial):")
for i, c in enumerate(coeficientes_monomial):
    print(f"a_{i} = {c:.6e}")

# Lagrange
def lagrange_basis(x, i, X):
    """Calcula el polinomio base de Lagrange L_i(x)."""
    L = 1.0
    n = len(X)
    for j in range(n):
        if j != i:
            L *= (x - X[j]) / (X[i] - X[j])
    return L

def lagrange_interpol(x, X, Y):
    """Evalúa el polinomio de Lagrange en el punto x."""
    n = len(X)
    P = 0.0
    for i in range(n):
        P += Y[i] * lagrange_basis(x, i, X)
    return P

def obtener_coeficientes_lagrange(T, C):
    n = len(T)
    coef_resultado = np.zeros(n)
    
    for i in range(n):
        coef_Li = np.array([1.0])
        for j in range(n):
            if j != i:
                denominador = T[i] - T[j]
                nuevo_factor = np.array([-T[j], 1.0]) / denominador
                coef_Li = np.convolve(coef_Li, nuevo_factor)
        coef_resultado += C[i] * coef_Li

    return coef_resultado

coef_lagrange = obtener_coeficientes_lagrange(T, C)

print("\nCoeficientes del polinomio interpolante (Lagrange, base monomial):")
for i, c in enumerate(coef_lagrange):
    print(f"a_{i} = {c:.6e}")

C_interp_lagrange = [lagrange_interpol(t, T, C) for t in T_suave]

# Gráfica
plt.figure(figsize=(10, 6))
plt.plot(T, C, 'ko', label="Datos originales")
plt.plot(T_suave, C_interp_newton, label="Interpolación Newton")
plt.plot(T_suave, C_interp_monomial, label="Interpolación Monomial")
plt.plot(T_suave, C_interp_lagrange, label="Interpolación Lagrange")
plt.xlabel("Temperatura (T)")
plt.ylabel("Calor específico (C)")
plt.title("Interpolación polinómica: Newton, Monomial y Lagrange")
plt.legend()
plt.grid(True)
plt.show()
