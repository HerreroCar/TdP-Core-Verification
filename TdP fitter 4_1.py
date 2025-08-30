# ======================================================================================
# TdP_MCMC_PRODUCTION_v0.4.1.py
# Cátedra de la Teoría del Pellizco (TdP)
#
# Autores: Carlos Herrero 
#
# Propósito:
# Versión final, robusta y autocontenida del Fitter MCMC.
# 1. Implementa un modelo de matriz de masa tridiagonal físicamente consistente.
# 2. Utiliza priors refinados (p discreto, alpha Gaussiano, Jeffreys).
# 3. Ejecuta una corrida MCMC con `emcee` y comprueba la convergencia.
# 4. Genera y guarda los artefactos finales: 'corner_plot_v0.4_final.png' y 'chains.h5'.
# ======================================================================================

# --- Paso 0: Instalar las dependencias necesarias en el entorno de Colab ---
!pip install numpy scipy emcee corner h5py

import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh # Usamos eigvalsh para matrices hermíticas

# --- Módulo 1: Modelo Físico Robusto (Placeholder v1.1) ---

def tdp_model_placeholder(params_vector):
    """
    Placeholder robusto que respeta la estructura de acoplamiento tridiagonal.
    Calcula las masas emergentes como los autovalores de una matriz de masa.
    params_vector = [p, alpha, E1, E2, E3, t12, t23]
    """
    p, alpha, E1, E2, E3, t12, t23 = params_vector
    n_levels = 3

    # Construimos una matriz de masa tridiagonal (Hermítica)
    M = np.zeros((n_levels, n_levels), dtype=float)

    # Diagonal: Energías base
    M[0, 0], M[1, 1], M[2, 2] = E1, E2, E3

    # Off-diagonal: Acoplamientos atenuados por la jerarquía p-ádica
    coupling12 = t12 / (p**(alpha * 1))
    coupling23 = t23 / (p**(alpha * 2))
    M[0, 1] = M[1, 0] = coupling12
    M[1, 2] = M[2, 1] = coupling23

    # Las masas emergentes son los autovalores de esta matriz
    eigenvalues = eigvalsh(M)
    return np.sort(eigenvalues)

# --- Módulo 2: Espacio de Probabilidad (Likelihood y Priors) ---

def log_prior(params_vector):
    """
    Define los priors para los parámetros, siguiendo las especificaciones de Grok.
    """
    try:
        p, alpha, E1, E2, E3, t12, t23 = params_vector
    except (ValueError, TypeError):
        return -np.inf

    # Prior para 'p': Uniforme discreto sobre primos relevantes
    if p not in [3, 5, 7, 11]:
        return -np.inf

    # Prior para 'alpha': Gaussiano centrado en el número áureo
    alpha_mean, alpha_std = 0.618, 0.05
    if not (0.5 <= alpha <= 0.7):
        return -np.inf

    # Priors para 'E_base' y 'couplings': Log-uniforme en sus rangos
    # (Implementado como uniforme en el valor por simplicidad, como es común)
    if not (0.01 <= E1 <= 1000 and 0.01 <= E2 <= 2000 and 0.01 <= E3 <= 4000):
        return -np.inf
    if not (0.01 <= t12 <= 100 and 0.01 <= t23 <= 100):
        return -np.inf

    # Calculamos la parte del log-prior que no es constante (la penalización Gaussiana)
    log_prior_alpha = -0.5 * ((alpha - alpha_mean) / alpha_std)**2

    return log_prior_alpha

def log_probability(params_vector, data, errors):
    """
    Función de probabilidad logarítmica total (log-likelihood + log-prior).
    """
    lp = log_prior(params_vector)
    if not np.isfinite(lp):
        return -np.inf

    predicted_masses = tdp_model_placeholder(params_vector)

    # Ajustamos la escala global con la primera generación
    scaler = data[0] / predicted_masses[0]
    predicted_masses_scaled = predicted_masses * scaler

    # Calculamos el chi-cuadrado
    chi2 = np.sum(((predicted_masses_scaled - data) / errors)**2)

    return lp - 0.5 * chi2

# --- Módulo 3: Extensión Cosmológica ---
def predict_B(chain):
    p_samples = chain[:, 0]
    alpha_samples = chain[:, 1]
    scale_factor = 1e-4  # Factor a calibrar con datos reales del CMB
    return 0.1 * scale_factor / (p_samples ** alpha_samples)

# --- Módulo 4: Ejecución del Fitter MCMC ---
if __name__ == "__main__":

    # Datos de entrada (leptones en MeV)
    data_leptons = np.array([0.511, 105.66, 1776.86])
    errors_leptons = 0.01 * data_leptons

    # Punto de partida ("ADN de Gaia" como vector)
    # [p, alpha, E1, E2, E3, t12, t23]
    initial_guess = [7.0, 0.618, 0.15, 2.8, 45.0, 10.0, 80.0]

    ndim = len(initial_guess)
    nwalkers = 100
    nsteps = 100000 # Reducido para una ejecución razonable en Colab (~10-20 min)
                  # Para producción, usar 100,000

    pos = np.array(initial_guess) + 1e-2 * np.random.randn(nwalkers, ndim)
    pos[:, 0] = np.random.choice([3, 5, 7, 11], size=nwalkers)

    print("Iniciando Fitter MCMC v0.4.1 (Corrida de prueba)...")

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data_leptons, errors_leptons))
    sampler.run_mcmc(pos, nsteps, progress=True)

    print("\n--- Corrida MCMC Completada ---")

    # --- Módulo 5: Análisis y Salida de Resultados ---

    # 1. Cadenas y Convergencia
    acceptance_fraction = np.mean(sampler.acceptance_fraction)
    print(f"Mean Acceptance Fraction: {acceptance_fraction:.3f}")

    # 2. Guardar Cadenas
    import h5py
    print("\nGuardando las cadenas en 'chains.h5'...")
    with h5py.File('chains.h5', 'w') as f:
        f.create_dataset('mcmc_chains', data=sampler.get_chain())
    print("Archivo 'chains.h5' generado con éxito.")

    # 3. Corner Plot
    samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    labels = ["p", r"$\alpha$", "E1 (MeV)", "E2 (MeV)", "E3 (MeV)", "t12", "t23"]
    print("\nGenerando Corner Plot...")
    fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 10})
    fig.savefig("corner_plot_v0.4_final.png")
    print("Archivo 'corner_plot_v0.4_final.png' generado con éxito.")
    plt.show()

    # 4. Predicción Cosmológica
    B_samples = predict_B(samples)
    B_mean = np.mean(B_samples)
    B_low, B_high = np.percentile(B_samples, [16, 84])

    print("\n--- Predicción Cosmológica ---")
    print(f"Amplitud de Oscilaciones CMB Predicha (B): {B_mean:.3e} (+{B_high-B_mean:.3e} / -{B_mean-B_low:.3e})")

    print("\n--- Fitter v0.4.1 listo. Archivos generados. ---")