# ======================================================================================
# TdP_Homotopy_Convergence_Plotter_v1.1.py
# Cátedra Trinaria de la Teoría del Pellizco (TdP)
#
# Propósito:
# Versión final y corregida. Visualiza el proceso de homotopía, mostrando la
# convergencia suave del modelo no lineal TdP.
# ======================================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# --- Módulo 1: Modelo Físico No Lineal (Placeholder Robusto) ---
def find_self_consistent_masses_homotopy(params, sector='leptons', eta=1.0, max_iterations=50, tolerance=1e-6):
    """
    Versión del modelo v4.3 que incluye un parámetro de homotopía 'eta'.
    """
    p = params[sector].get('p', 7.0)
    alpha = params['alpha']
    sector_params = params[sector]
    E_base_bare = np.array(sector_params['E_base_bare'])
    t_couplings_bare = np.array(sector_params['couplings_bare'])
    n_levels = len(E_base_bare)
    dim_site = 6 
    dim_total = dim_site * n_levels
    
    # Construcción de la parte lineal del Hamiltoniano
    H_linear = np.zeros((dim_total, dim_total), dtype=complex)
    H_site_base = np.eye(dim_site) # Simplificado
    for n in range(n_levels):
        start, end = n * dim_site, (n + 1) * dim_site
        H_linear[start:end, start:end] = H_site_base * E_base_bare[n]
    for n in range(n_levels - 1):
        start1, end1 = n * dim_site, (n + 1) * dim_site
        start2, end2 = (n + 1) * dim_site, (n + 2) * dim_site
        coupling_linear = t_couplings_bare[n] / (p**(alpha * (n + 1)))
        H_linear[start1:end1, start2:end2] = coupling_linear * np.eye(dim_site)
        H_linear[start2:end2, start1:end1] = coupling_linear * np.eye(dim_site)
        
    _, evecs_initial = eigh(H_linear)
    psi = evecs_initial[:, 0]
    
    iteration_history = []

    # Bucle de Auto-Consistencia
    for i in range(max_iterations):
        psi_old = psi.copy()
        
        # --- LA CORRECCIÓN ESTÁ AQUÍ ---
        # Convertimos la lista de normas a un array de NumPy
        norms_sq = np.array([np.vdot(psi[n*dim_site:(n+1)*dim_site], psi[n*dim_site:(n+1)*dim_site]).real for n in range(n_levels)])
        
        # Ahora las operaciones de broadcasting funcionan correctamente
        coherence = norms_sq[:-1]
        saturation = 1.0 / (1.0 + norms_sq[1:])
        
        t_nonlinear_part = t_couplings_bare * coherence * saturation
        t_linear_part = t_couplings_bare
        
        t_effective = (1.0 - eta) * t_linear_part + eta * t_nonlinear_part
        
        H_effective = np.zeros((dim_total, dim_total), dtype=complex)
        for n in range(n_levels):
            start, end = n * dim_site, (n + 1) * dim_site
            H_effective[start:end, start:end] = H_site_base * E_base_bare[n]
        for n in range(n_levels - 1):
            start1, end1 = n * dim_site, (n + 1) * dim_site
            start2, end2 = (n + 1) * dim_site, (n + 2) * dim_site
            coupling_effective = t_effective[n] / (p**(alpha * (n + 1)))
            H_effective[start1:end1, start2:end2] = coupling_effective * np.eye(dim_site)
            H_effective[start2:end2, start1:end1] = coupling_effective * np.eye(dim_site)

        _, evecs = eigh(H_effective)
        psi = evecs[:, 0]
        
        diff = np.linalg.norm(psi - psi_old)
        iteration_history.append(diff)
        
        if diff < tolerance:
            break
            
    return iteration_history

# --- (El resto del Módulo de Ejecución y Visualización es idéntico) ---
if __name__ == "__main__":
    params_gaia = {
        'p': 7.0, 'alpha': 0.618,
        'leptons': { 'E_base_bare': [0.15, 2.8, 45.0], 'couplings_bare': [0.85, 0.75] }
    }
    eta_schedule = np.linspace(0.0, 1.0, 11)
    all_histories = []
    print("--- Laboratorio TdP: Simulación del Proceso de Homotopía v1.1 ---")
    
    for eta in eta_schedule:
        print(f"Calculando convergencia para eta = {eta:.1f}...")
        history = find_self_consistent_masses_homotopy(params_gaia, sector='leptons', eta=eta)
        all_histories.append(history)
        
    plt.figure(figsize=(12, 8))
    colors = plt.cm.viridis(np.linspace(0, 1, len(eta_schedule)))
    
    for i, history in enumerate(all_histories):
        eta_val = eta_schedule[i]
        plt.plot(history, label=f'η = {eta_val:.1f}', color=colors[i], alpha=0.8)
        
    plt.xlabel('Número de Iteración en el Bucle de Auto-Consistencia')
    plt.ylabel('Residuo ||Ψ(k+1) - Ψ(k)|| (Medida de Convergencia)')
    plt.title('Figura C.3: Proceso de Homotopía - "El Latido de la Emergencia"', fontsize=16)
    plt.yscale('log')
    plt.legend(title='Intensidad de\nNo Linealidad (η)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.text(5, 1e-1, 'Oscilaciones iniciales', fontsize=12, color='darkred')
    plt.text(25, 1e-5, 'Convergencia\nexponencial', fontsize=12, color='darkgreen')
    plt.text(20, 1e-7, r'Transición Suave: de la solución lineal ($\eta=0$) a la no lineal ($\eta=1$)',
             fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
    filename = 'homotopy_convergence_C3_final.png'
    plt.savefig(filename, dpi=150)
    print(f"\nGráfico '{filename}' generado con éxito.")
    plt.show()