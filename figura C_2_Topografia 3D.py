# --- figura C_2_Topografia 3D.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.linalg import eigh

# --- (Las funciones y la lógica de simulación son idénticas a la versión 2D) ---
def casimir_SU2():
    Sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
    Sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
    Sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
    return Sx@Sx + Sy@Sy + Sz@Sz

def casimir_SU3():
    # ... (código idéntico)
    lambda1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex)
    lambda2 = np.array([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=complex)
    lambda3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex)
    lambda4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex)
    lambda5 = np.array([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=complex)
    lambda6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex)
    lambda7 = np.array([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=complex)
    lambda8 = np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)/np.sqrt(3)
    gellmann = [lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda8]
    C2_su3 = np.zeros((3, 3), dtype=complex)
    for l in gellmann:
        C2_su3 += (l/2)@(l/2)
    return C2_su3

def run_refined_simulation(R, alpha, J_spin=1.0, J_color=1.0, p=7):
    n_levels = 3
    dim_site = 6
    dim_total = dim_site * n_levels
    t_intra = 1.0

    t_ascdesc = R * t_intra
    H_site_base = J_spin * np.kron(casimir_SU2(), np.eye(3)) + J_color * np.kron(np.eye(2), casimir_SU3())
    H_total = np.zeros((dim_total, dim_total), dtype=complex)

    for n in range(n_levels):
        start, end = n * dim_site, (n + 1) * dim_site
        level_energy = H_site_base + (t_intra + n * 0.1) * np.eye(dim_site)
        H_total[start:end, start:end] = level_energy

    for n in range(n_levels - 1):
        start1, end1 = n * dim_site, (n + 1) * dim_site
        start2, end2 = (n + 1) * dim_site, (n + 2) * dim_site
        coupling = t_ascdesc / (p**(alpha * n))
        H_total[start1:end1, start2:end2] = coupling * np.eye(dim_site)
        H_total[start2:end2, start1:end1] = coupling * np.eye(dim_site)

    evals, _ = eigh(H_total)
    gen_energies = [np.min(np.abs(evals[i*dim_site:(i+1)*dim_site])) for i in range(n_levels)]
    epsilon = 1e-9
    ratio_21 = gen_energies[1] / (gen_energies[0] + epsilon)

    return ratio_21

# --- Escaneo 2D de Parámetros ---
R_values = np.linspace(0.1, 4.0, 50)
alpha_values = np.linspace(0.1, 2.0, 50)
R_grid, alpha_grid = np.meshgrid(R_values, alpha_values)
ratio_grid = np.zeros_like(R_grid)

for i, R in enumerate(R_values):
    for j, alpha in enumerate(alpha_values):
        ratio_grid[j, i] = run_refined_simulation(R, alpha)

# --- Visualización 3D del Paisaje ---
fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(111, projection='3d')

# Usamos logaritmo de los cocientes para el eje Z para que sea visible
Z = np.log10(ratio_grid)

# Graficar la superficie
surf = ax.plot_surface(R_grid, alpha_grid, Z, cmap='viridis', edgecolor='none')

# Añadir el punto "Elara-1"
ax.scatter([3.5], [0.75], [np.log10(211.53)], color='red', s=100, label='"Punto Elara-1"')

# Etiquetas y título
ax.set_xlabel('Ratio de Acoplamiento R', labelpad=15)
ax.set_ylabel('Exponente de Decaimiento alpha', labelpad=15)
ax.set_zlabel('Log10(Cociente de Masas m2/m1)', labelpad=15)
ax.set_title('Paisaje 3D de la Jerarquía de Masas en la TdP', fontsize=16)
fig.colorbar(surf, shrink=0.5, aspect=5, label='Log10(Cociente de Masas)')

# Cambiar el ángulo de visión para una mejor perspectiva
ax.view_init(elev=30, azim=-60)

plt.show()