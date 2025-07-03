import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parámetros del modelo
u = 10
n = 0.01
beta = 0.04  # Supuesto beta en estado endémico

# Equilibrio endémico (estimado)
Xs_eq = u / (beta + n)
Xi_eq = 80
Xr_eq = 60

# Sistema reducido a 3 variables: Xs, Xi, Xr
def tuberculosis_deriv(X, t):
    Xs, Xi, Xr = X
    dXs = u - (beta + n) * Xs
    dXi = beta * Xs - n * Xi - 0.05 * Xi
    dXr = 0.05 * Xi - n * Xr
    return np.array([dXs, dXi, dXr])

# Tiempo para integración
t = np.linspace(0, 200, 500)

# Generar múltiples condiciones iniciales alrededor del equilibrio
initial_conditions = [
    [Xs_eq * 0.9, Xi_eq * 0.8, Xr_eq * 0.85],
    [Xs_eq * 1.1, Xi_eq * 1.2, Xr_eq * 1.15],
    [Xs_eq * 1.0, Xi_eq * 0.6, Xr_eq * 1.1],
    [Xs_eq * 0.8, Xi_eq * 1.4, Xr_eq * 0.9],
    [Xs_eq * 1.05, Xi_eq * 1.0, Xr_eq * 1.3],
    [Xs_eq * 1.2, Xi_eq * 0.9, Xr_eq * 1.0]
]

# Soluciones de las trayectorias
trajectories = []
for X0 in initial_conditions:
    X = np.zeros((len(t), 3))
    X[0] = X0
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        X[i] = X[i-1] + tuberculosis_deriv(X[i-1], t[i-1]) * dt
    trajectories.append(X)

# Graficar
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

colors = ['r', 'g', 'b', 'm', 'c', 'orange']
for X, color in zip(trajectories, colors):
    ax.plot(X[:, 0], X[:, 1], X[:, 2], label=f"CI: [{X[0,0]:.1f}, {X[0,1]:.1f}, {X[0,2]:.1f}]", color=color)

ax.set_xlabel('Susceptibles (Xs)')
ax.set_ylabel('Infectados (Xi)')
ax.set_zlabel('Recuperados (Xr)')
ax.set_title('Trayectorias hacia el Equilibrio Endémico (EE)')
ax.legend()
plt.tight_layout()
plt.show()
