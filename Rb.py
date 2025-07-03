import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Parámetros fijos
n = 0.01477
c = 0.85
g = 0.05
h = 0.98
gamma = 0.3
r = 0.16
d1 = 0.02
d2 = 0.006
m1 = 0.5
m2 = 0.3
p = 0.7
k = 0.4
b = 2.001
a1 = 0.3
a2 = 0.3

# Función para calcular Rb
def calcular_rb(a1, a2, b, c, g, k, h, gamma, r, p, m1, m2, n, d1, d2):
    z1 = n + d1 + h
    z2 = n + gamma
    z3 = n + d2 + r
    z4 = n + gamma
    num = (1 - a1) * b * (c * g + (1 - c) * z1) * ((1 - a2) * z4 * (z3 + m1 * k * h) + (1 - k) * h * m2 * z3)
    den = z3 * z1 * z2 * z4 - (1 - p) * r * g * (1 - k) * h
    return num / den if den != 0 else np.nan

# Resolución de grilla
res = 40

# Variables de las 4 gráficas
a1_vals = np.linspace(0, 1, res)
a2_vals = np.linspace(0, 1, res)
gamma_vals = np.linspace(0.01, 1, res)
b_vals = np.linspace(0.5, 3, res)
p_vals = np.linspace(0, 1, res)
k_vals = np.linspace(0, 1, res)
g_vals = np.linspace(0.01, 0.2, res)

# Superficies de Rb
rb_a1_a2 = np.array([[calcular_rb(a1, a2, b, c, g, k, h, gamma, r, p, m1, m2, n, d1, d2) for a2 in a2_vals] for a1 in a1_vals])
rb_gamma_b = np.array([[calcular_rb(a1, a2, b_, c, g, k, h, gamma_, r, p, m1, m2, n, d1, d2) for b_ in b_vals] for gamma_ in gamma_vals])
rb_p_k = np.array([[calcular_rb(a1, a2, b, c, g, k_, h, gamma, r, p_, m1, m2, n, d1, d2) for k_ in k_vals] for p_ in p_vals])
rb_g_k = np.array([[calcular_rb(a1, a2, b, c, g_, k_, h, gamma, r, p, m1, m2, n, d1, d2) for k_ in k_vals] for g_ in g_vals])

# Crear subplots 2x2
fig = make_subplots(
    rows=2, cols=2,
    specs=[[{'type': 'surface'}, {'type': 'surface'}],
           [{'type': 'surface'}, {'type': 'surface'}]],
    subplot_titles=[
        "Rb vs a1 y a2",
        "Rb vs γ y b",
        "Rb vs p y k",
        "Rb vs g y k"
    ]
)

# Añadir las 4 superficies
fig.add_trace(go.Surface(z=rb_a1_a2, x=a1_vals, y=a2_vals, colorscale='Viridis', showscale=False),
              row=1, col=1)

fig.add_trace(go.Surface(z=rb_gamma_b, x=gamma_vals, y=b_vals, colorscale='Cividis', showscale=False),
              row=1, col=2)

fig.add_trace(go.Surface(z=rb_p_k, x=p_vals, y=k_vals, colorscale='Plasma', showscale=False),
              row=2, col=1)

fig.add_trace(go.Surface(z=rb_g_k, x=g_vals, y=k_vals, colorscale='Inferno', showscale=False),
              row=2, col=2)

# Títulos y ejes
fig.update_layout(
    title="Comportamiento tridimensional del número umbral Rb",
    height=900,
    width=1200,
    margin=dict(l=10, r=10, t=40, b=10)
)

# Ejes individuales (opcional: puedes quitar para limpieza)
fig.update_scenes(
    xaxis_title="a1", yaxis_title="a2", zaxis_title="Rb", row=1, col=1
)
fig.update_scenes(
    xaxis_title="γ", yaxis_title="b", zaxis_title="Rb", row=1, col=2
)
fig.update_scenes(
    xaxis_title="p", yaxis_title="k", zaxis_title="Rb", row=2, col=1
)
fig.update_scenes(
    xaxis_title="g", yaxis_title="k", zaxis_title="Rb", row=2, col=2
)

# Mostrar
fig.show()