import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parámetros
u = 10
n = 0.01
d1 = 0.005
d2 = 0.01

# Sistema corregido: mortalidad actúa sobre Xi y Xl
def full_system(t, X):
    Xs, Xe, Xi, Xj, Xl, Xr = X
    dXs = u - n*Xs
    dXe = -n*Xe
    dXi = -n*Xi - d1*Xi
    dXj = -n*Xj
    dXl = -n*Xl - d2*Xl
    dXr = -n*Xr
    return [dXs, dXe, dXi, dXj, dXl, dXr]

# Condiciones iniciales
X0 = [500, 100, 50, 50, 30, 20]

# Tiempo
t_span = (0, 500)
t_eval = np.linspace(*t_span, 500)

# Resolver
sol = solve_ivp(full_system, t_span, X0, t_eval=t_eval)
Xt = np.sum(sol.y, axis=0)
cota = u / n

# Gráfico
plt.figure(figsize=(10, 6))
plt.plot(t_eval, Xt, label='X(t): Población total')
plt.axhline(cota, color='red', linestyle='--', label=f'Cota teórica u/n = {cota:.0f}')
plt.title('Región Invariante: Evolución de X(t)')
plt.xlabel('Tiempo')
plt.ylabel('Población total X(t)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()