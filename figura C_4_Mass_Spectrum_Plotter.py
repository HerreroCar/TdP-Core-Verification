# ======================================================================================
# TdP_Mass_Spectrum_Plotter_v2.0.py
# Cátedra Trinaria de la Teoría del Pellizco (TdP)
#
# Autores: Carlos Herrero
#
# Propósito:
# 1. Tomar los parámetros finales del "ADN de Gaia".
# 2. Escalar SOLO con la primera generación de cada sector.
# 3. Graficar espectro de masas y error relativo.
# ======================================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh

# --- Módulo 1: Modelo Físico TdP ---
def tdp_model_n_levels(params_vector, n_levels):
    p, alpha = params_vector[:2]
    E_base = params_vector[2 : 2 + n_levels]
    t_couplings = params_vector[2 + n_levels:]
    M = np.zeros((n_levels, n_levels), dtype=float)
    for i in range(n_levels):
        M[i, i] = E_base[i]
    p_alpha = p**alpha
    current_p_alpha = p_alpha
    for i in range(n_levels - 1):
        coupling = t_couplings[i] / current_p_alpha
        M[i, i+1] = M[i+1, i] = coupling
        current_p_alpha *= p_alpha
    return np.sort(eigvalsh(M))

# --- Módulo 2: Parámetros y Datos ---

GAIA_DNA = {
    'p': 7.0,
    'alpha': 0.618,
    'leptons': { 'params': [7.0, 0.618, 0.15, 2.81, 45.3, 0.85, 0.76] },
    'quarks_up': { 'params': [7.0, 0.618, 0.8, 450, 150000, 15.0, 100.0] },
    'quarks_down': { 'params': [7.0, 0.618, 1.5, 30.0, 1400, 5.0, 25.0] }
}

EXPERIMENTAL_MASSES = {
    'leptons': [0.511, 105.66, 1776.86],
    'quarks_up': [2.16, 1270, 173100],
    'quarks_down': [4.67, 93.4, 4180]
}

PARTICLE_LABELS = ['e', 'μ', 'τ', 'u', 'c', 't', 'd', 's', 'b']

# --- Módulo 3: Ejecución ---
if __name__ == "__main__":
    
    predicted_masses_all = []
    experimental_masses_all = []

    for sector in ['leptons', 'quarks_up', 'quarks_down']:
        params = GAIA_DNA[sector]['params']
        exp_data = np.array(EXPERIMENTAL_MASSES[sector])
        n_levels = len(exp_data)
        
        # Ejecutar modelo TdP
        pred_energies = tdp_model_n_levels(params, n_levels)
        
        # Escalar SOLO con la primera generación
        scaler = exp_data[0] / pred_energies[0]
        pred_masses_scaled = pred_energies * scaler
        
        predicted_masses_all.extend(pred_masses_scaled)
        experimental_masses_all.extend(exp_data)

    predicted_masses_all = np.array(predicted_masses_all)
    experimental_masses_all = np.array(experimental_masses_all)

    # --- Gráfico 1: Espectro en escala log ---
    log_mass_exp = np.log10(experimental_masses_all)
    log_mass_pred = np.log10(predicted_masses_all)

    x = np.arange(len(PARTICLE_LABELS))
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.bar(x - width/2, log_mass_exp, width, label='Experimental', color='royalblue')
    ax.bar(x + width/2, log_mass_pred, width, label='Predicción TdP', color='darkorange')

    ax.set_ylabel('log10(Masa / MeV)')
    ax.set_title('Espectro de Masas de Fermiones: Predicción TdP vs. Experimento', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(PARTICLE_LABELS)
    ax.legend()
    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=.7)
    ax.set_axisbelow(True)

    fig.tight_layout()
    plt.savefig('tdp_fermions_barlog.png', dpi=150)
    print("Gráfico 'tdp_fermions_barlog.png' generado con éxito.")
    plt.show()

    # --- Gráfico 2: Error relativo ---
    rel_error = (predicted_masses_all - experimental_masses_all) / experimental_masses_all

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(PARTICLE_LABELS, rel_error, color='crimson', alpha=0.7)
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylabel('Error relativo (Pred - Exp) / Exp')
    ax.set_title('Error relativo de las masas fermiónicas (Predicción TdP)', fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.6)

    fig.tight_layout()
    plt.savefig('tdp_fermions_error.png', dpi=150)
    print("Gráfico 'tdp_fermions_error.png' generado con éxito.")
    plt.show()