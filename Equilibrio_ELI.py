import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
b = 0.3
a1 = 0.2
a2 = 0.2
m1 = 0.5
m2 = 0.4
gamma = 0.05

# Tasa de infección
def beta(X_i, X_j, X_l, X_total):
    return (b * (1 - a1) / X_total) * ((1 - a2) * X_i + m1 * (1 - a2) * X_j + m2 * X_l)

# Sistema de EDOs
def tuberculosis_model(t, X):
    X_s, X_e, X_i, X_j, X_l, X_r = X
    X_total = np.sum(X)
    B = beta(X_i, X_j, X_l, X_total)

    dXs = u - (B + n) * X_s
    dXe = c * B * X_s - (n + g) * X_e + (1 - p) * r * X_l
    dXi = (1 - c) * B * X_s + g * X_e - (n + d1 + h) * X_i
    dXj = k * h * X_i - (n + y) * X_j
    dXl = (1 - k) * h * X_i - (n + d2 + r) * X_l
    dXr = y * X_j + p * r * X_l - n * X_r

    return [dXs, dXe, dXi, dXj, dXl, dXr]

# Condiciones iniciales y tiempo
t_eval = np.linspace(0, 200, 1000)
initial_conditions = [
    [800, 0, 0, 0, 0, 0],
    [400, 10, 5, 3, 2, 1],
    [600, 20, 30, 10, 5, 0],
    [300, 0, 20, 15, 8, 0],
]

# Punto de equilibrio libre de infección
X_eli = [u / n, 0, 0, 0, 0, 0]

# Integración y gráficas
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

for ic in initial_conditions:
    sol = solve_ivp(tuberculosis_model, [0, 200], ic, t_eval=t_eval)
    X_s, X_e, X_i = sol.y[0], sol.y[1], sol.y[2]
    ax.plot3D(X_s, X_e, X_i, label=f"IC: {ic}")

# Punto ELI
ax.scatter(*X_eli[:3], color='blue', s=60, marker='o', label='Equilibrio Libre (ELI)', zorder=5)


# Etiquetas y estilo
ax.set_xlabel("$X_s$ (Susceptibles)")
ax.set_ylabel("$X_e$ (Expuestos)")
ax.set_zlabel("$X_i$ (Infectados)")
ax.set_title("Trayectorias hacia los Puntos de Equilibrio")
ax.legend()
plt.tight_layout()
plt.show()
