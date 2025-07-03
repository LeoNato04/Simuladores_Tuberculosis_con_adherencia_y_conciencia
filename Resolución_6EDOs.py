import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

# Crear carpeta de salida
output_dir = "/mnt/data/simulacion_tuberculosis_final_variaciones"
os.makedirs(output_dir, exist_ok=True)

# Parámetros base
params_base = {
    'u': 10, 'n': 0.01, 'b': 0.5, 'a1': 0.1, 'a2': 0.1,
    'c': 0.3, 'g': 0.05, 'h': 0.1, 'k': 0.6, 'gamma': 0.1,
    'r': 0.05, 'p': 0.7, 'd1': 0.005, 'd2': 0.01, 'm1': 0.5, 'm2': 0.3
}

# Condiciones iniciales
Y0 = [710800, 8000, 217, 510, 768, 0]
T = 200
dt = 1
t_eval = list(range(0, T + 1, dt))

# Sistema de EDOs
def sistema(t, y, p):
    S, E, I, J, L, R = y
    X = S + E + I + J + L + R
    beta = p['b'] * (1 - p['a1']) * ((1 - p['a2']) * I + p['m1'] * (1 - p['a2']) * J + p['m2'] * L) / X
    dS = p['u'] - (beta + p['n']) * S
    dE = p['c'] * beta * S - (p['n'] + p['g']) * E + (1 - p['p']) * p['r'] * L
    dI = (1 - p['c']) * beta * S + p['g'] * E - (p['n'] + p['d1'] + p['h']) * I
    dJ = p['k'] * p['h'] * I - (p['n'] + p['gamma']) * J
    dL = (1 - p['k']) * p['h'] * I - (p['n'] + p['d2'] + p['r']) * L
    dR = p['gamma'] * J + p['p'] * p['r'] * L - p['n'] * R
    return [dS, dE, dI, dJ, dL, dR]

# Variables del modelo
variables = ['S', 'E', 'I', 'J', 'L', 'R']
edo_parametros = {
    'S': ['b', 'a1', 'a2', 'n'],
    'E': ['c', 'g', 'p'],
    'I': ['c', 'g', 'd1', 'h'],
    'J': ['k', 'h', 'gamma'],
    'L': ['k', 'r', 'd2'],
    'R': ['gamma', 'p', 'r']
}
variaciones = [0.8, 1.2]

# Simulación con variaciones por EDO
for i, var in enumerate(variables):
    plt.figure(figsize=(10, 6))

    # Solución base
    sol_base = solve_ivp(lambda t, y: sistema(t, y, params_base), [0, T], Y0, t_eval=t_eval)
    plt.plot(sol_base.t, sol_base.y[i], label=f'{var} base', color='black', linestyle='--')

    for param in edo_parametros[var]:
        for factor in variaciones:
            mod_params = params_base.copy()
            mod_params[param] = round(mod_params[param] * factor, 5)
            sol_mod = solve_ivp(lambda t, y: sistema(t, y, mod_params), [0, T], Y0, t_eval=t_eval)
            label = f"{param} = {mod_params[param]}"
            plt.plot(sol_mod.t, sol_mod.y[i], label=label)

    plt.title(f'Variaciones paramétricas en {var}(t)')
    plt.xlabel('Tiempo (días)')
    plt.ylabel(f'{var}(t)')
    plt.legend(fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{var}_variaciones.png"))
    plt.close()

output_dir
