import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# Parámetros del modelo
u = 100
n = 0.1
g = 0.05
d1 = 0.03
d2 = 0.02
r = 0.04
y = 0.05
h = 0.08
k = 0.6
c = 0.7
p = 0.8
b = 0.9 
a1 = 0.2
a2 = 0.2
m1 = 0.5
m2 = 0.4

# Cálculo de variables intermedias
z1 = n + g
z2 = n + d1 + h
z3 = n + y      # y ≡ gamma
z4 = n + d2 + r

# Numerador y denominador del Rb
numerador = (1 - a1) * b * (c * g + (1 - c) * z1)
factor_infeccioso = (1 - a2) * z4 * (z3 + m1 * k * h) + (1 - k) * h * m2 * z3
denominador = z3 * (z1 * z2 * z4 - (1 - p) * r * (1 - k) * h)

# Cálculo final de Rb
Rb = (numerador * factor_infeccioso) / denominador

print(f"El valor de Rb con los parámetros actuales es: {Rb:.4f}")

# Tasa de infección beta
def beta(X_i, X_j, X_l, X_total):
    return (b * (1 - a1) / X_total) * ((1 - a2) * X_i + m1 * (1 - a2) * X_j + m2 * X_l)

# Sistema reducido a Xe y Xi (manteniendo valores bajos en otros compartimentos)
def tuberculosis_plane(t, X):
    X_e, X_i = X
    X_s = u / n  # equilibrio libre de infección
    X_j = 0.1
    X_l = 0.1
    X_r = 0.1
    X_total = X_s + X_e + X_i + X_j + X_l + X_r
    B = beta(X_i, X_j, X_l, X_total)
    
    dXe = c * B * X_s - (n + g) * X_e + (1 - p) * r * X_l
    dXi = (1 - c) * B * X_s + g * X_e - (n + d1 + h) * X_i
    return [dXe, dXi]

# Crear malla
x_vals = np.linspace(-1, 5, 25)
y_vals = np.linspace(-1, 5, 25)
Xe, Xi = np.meshgrid(x_vals, y_vals)

# Campo vectorial
dXe_vals = np.zeros_like(Xe)
dXi_vals = np.zeros_like(Xi)

for i in range(Xe.shape[0]):
    for j in range(Xe.shape[1]):
        dx = tuberculosis_plane(0, [Xe[i, j], Xi[i, j]])
        dXe_vals[i, j] = dx[0]
        dXi_vals[i, j] = dx[1]

# Graficar campo vectorial y líneas de flujo
plt.figure(figsize=(8, 6))
plt.streamplot(Xe, Xi, dXe_vals, dXi_vals, color='crimson', density=1.2)
plt.quiver(Xe, Xi, dXe_vals, dXi_vals, color='black', alpha=0.4)

plt.xlabel("$X_e$ (Expuestos)")
plt.ylabel("$X_i$ (Infectados)")
plt.title("Campo Vectorial y Líneas de Flujo (Plano $X_e$-$X_i$)")
plt.grid(True)
plt.axhline(0, color='blue', lw=0.5)
plt.axvline(0, color='blue', lw=0.5)
plt.tight_layout()
plt.show()
