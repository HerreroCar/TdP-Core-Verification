# --- figura C_1_Topografia 2D.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.linalg import eigh # ¡LA CORRECCIÓN ESTÁ AQUÍ! Importamos 'eigh' explícitamente.

# --- (Las funciones de Casimir son las mismas) ---
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

# --- Modelo Refinado ---
def run_refined_simulation(R, alpha, J_spin=1.0, J_color=1.0, p=7):
    n_levels = 3
    dim_site = 6
    dim_total = dim_site * n_levels
    t_intra = 1.0 # Fijamos la escala de energía base

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

# --- Visualización del Paisaje ---
plt.figure(figsize=(12, 9))
contour = plt.contourf(R_grid, alpha_grid, ratio_grid, levels=np.logspace(0, 3, 20), cmap='viridis', norm=LogNorm())
cbar = plt.colorbar(contour)
cbar.set_label('Cociente de Masas Emergentes $m_2/m_1$', fontsize=12)
plt.axhline(y=0.75, color='red', linestyle='--', linewidth=2)
plt.axvline(x=3.5, color='red', linestyle='--', linewidth=2)
plt.text(3.55, 0.77, '"Punto Elara-1"', color='red', fontsize=12, bbox=dict(facecolor='white', alpha=0.7))
plt.xlabel(r'Ratio de Acoplamiento $R = t_{ascdesc} / t_{intra}$', fontsize=14)
plt.ylabel(r'Exponente de Decaimiento $\alpha$', fontsize=14)
plt.title('Paisaje de la Jerarquía de Masas en la TdP', fontsize=16)
plt.show()