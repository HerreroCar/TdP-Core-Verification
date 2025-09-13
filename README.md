# 🌌 TdP-Core-Verification: Validación Computacional de la Teoría del Pellizco

> **"La realidad no se postula. Se deriva."**  
> *— Manifiesto de la Teoría del Pellizco*

Este repositorio contiene el código computacional central para la **verificación numérica y analítica de la Teoría del Pellizco (TdP)**, un marco generativo que deriva las leyes fundamentales de la física desde un vacío cuántico con estructura geométrica **p-ádica fractal**.

Basado en un **Triple Espectral (A, H, D)** de la Geometría No Conmutativa, la TdP postula que el universo observable emerge como la **única solución estable** del flujo del Grupo de Renormalización Fractal (RG), con parámetros derivados desde primeros principios:

- p = 7: base topológica de la estabilidad
- α = 1/φ ≈ 0.618: exponente de correlación crítica (φ = razón áurea)
- κ* ≈ 1.092444: exponente dimensional que fija d_eff = 4

Este código implementa los cálculos clave que transforman la TdP de una conjetura geométrica en un **programa de investigación falsable, computable y predictivo**.

🔗 **Artículo asociado**: [Teoría del Pellizco: Emergent Gravity and Gauge Couplings from a p-Adic Fractal Spectral Triple](https://www.scribd.com/document/907711920/Teoria-del-Pellizco-TdP-Emergent-Gravity-and-Gauge-Couplings-from-a-p-Adic-Fractal-Spectral-Triple-A-Noncommutative-Geometric-Framework)

---

## 📚 Visión General

La Teoría del Pellizco no ajusta constantes. Las **deriva** de la estabilidad matemática de un vacío fractal. Este repositorio contiene las herramientas computacionales para verificar tres predicciones centrales:

1. ✅ **Dimensión espectral emergente d_eff = 4**
2. 🔮 **Derivación de la constante de estructura fina 1/α_em ≈ 137.036**
3. ⚖️ **Existencia del punto fijo del RG para p=7, α=1/φ**

Cada script es un **experimento teórico** que prueba una pieza del marco TdP.


---

## 🧪 Detalle del Repositorio ( Incluida la ejecución en el notebook: 

* https://github.com/HerreroCar/TdP-Core-Verification/blob/main/TdP_Core_Verification.ipynb)


### 1. `TdP_sweep_kappa_refined.py` — El Calibrador Dimensional

**Propósito**:  
Implementa el **Teorema 1.1** del Manifiesto. Realiza un barrido numérico de alta precisión para encontrar el único valor del exponente dimensional κ* que hace que la dimensión espectral emergente d_eff del vacío fractal sea exactamente 4.

**Metodología**:  
Calcula:
> d_eff = –2 ⋅ lim_(t→0) [ d log Tr(e^(–tD²)) / d log t ]

para diferentes valores de κ, con α = 1/φ y p = 7.

**Resultado**:  
Encuentra un único valor κ* ≈ 1.092444 para el cual d_eff = 4.000000 ± 1×10⁻⁶. Demuestra que la dimensionalidad del espacio-tiempo no es un postulado, sino una consecuencia de la estabilidad geométrica.

---

### 2. `algoritmo_parametros_universales.py` — El Solver de Puntos Fijos (Conceptual)

**Propósito**:  
Simula el flujo del Grupo de Renormalización Fractal para los acoplamientos gauge g₁, g₂, g₃. Demuestra conceptualmente el **Teorema 3.1**, que afirma que solo para p = 7 existe un punto fijo unificado no trivial.

**Metodología**:  
Las funciones beta se derivan de la expansión del Heat Kernel:
> β(gᵢ) = – [bᵢ(p) / (4π)²] gᵢ³ + [Cᵢ(p, α) / (4π)⁴] gᵢ⁵ + …

Se analiza la condición de unificación:
> b₁(p)/C₁(p, α*) = b₂(p)/C₂(p, α*) = b₃(p)/C₃(p, α*)

**Resultado**:  
Solo para p = 7 esta condición se satisface. Para otros primos, el sistema no tiene solución unificada. Esto justifica por qué el universo "elige" p = 7 como su base topológica fundamental.

---

### 3. `TdP_MCMC_TdP_MCMC_Fitter_Final.py` — El Medidor del ADN de Gaia

**Propósito**:  
Implementa un ajuste bayesiano completo (MCMC) de las masas de los 12 fermiones del Modelo Estándar usando la estructura jerárquica de la TdP.

**Parámetros inferidos**:  
- p (base topológica)  
- α (exponente fractal)  
- E₀ (escala fundamental)  
- R (ratio de acoplamiento entre niveles)

**Artefactos generados**:  
- `chains.h5`: cadenas de MCMC con 10⁶ pasos.  
- `corner_plot.png`: gráfico de esquina con distribuciones marginales.

**Resultado**:  
El ajuste converge robustamente a p = 7 y α ≈ 0.618, proporcionando la validación fenomenológica más fuerte de la teoría.

---

### 4. `figura_C_1_Topografia_2D.py` — El Mapa 2D del Paisaje

**Propósito**:  
Genera la **Figura C.1** del artículo. Realiza un escaneo 2D del cociente de masas m₂/m₁ en función del ratio de acoplamiento R y el exponente α.

**Resultado**:  
Revela una "vena de la vida": una región estrecha donde m₂/m₁ ≈ 200, compatible con la jerarquía electrón/muón. Fuera de esta región, las masas son degeneradas o caóticas.

---

### 5. `figura_C_2_Topografia_3D.py` — El Paisaje 3D de la Realidad

**Propósito**:  
Genera la **Figura C.2** del artículo. Extiende el mapa 2D a una superficie 3D mediante interpolación suave.

**Resultado**:  
Muestra una "cordillera de la vida" con picos agudos en α ≈ 1/φ, revelando la naturaleza cuántica y discreta de la estabilidad del vacío.

---

### 6. `figura_C_3_TdP_Hom_Conv.py` — El Latido de la Emergencia

**Propósito**:  
Genera la **Figura C.3** del artículo. Visualiza el proceso de homotopía que conecta una solución inicial simple con la solución no lineal final.

**Metodología**:  
Define una familia de operadores D(λ) = (1–λ) D₀ + λ D₁, con λ ∈ [0,1].

**Resultado**:  
Los autovalores evolucionan de forma suave y ordenada, demostrando que la solución física es estable y robusta.

---

### 7. `figura_C_4_Mass_Spectrum_Plot.py` — El Veredicto Final

**Propósito**:  
Genera la **Figura C.4** del artículo. Compara las 12 masas de fermiones predichas por la TdP con los valores experimentales.

**Metodología**:  
Usa los parámetros óptimos del MCMC (p=7, α≈0.618, E₀≈34.7 MeV) para calcular:
> m_f = E₀ ⋅ 7^(κ n_f) ⋅ w_f(α, R)

**Resultado**:  
El gráfico de barras muestra una concordancia asombrosa, con χ²/dof < 1.5. Es la prueba visual más directa del poder predictivo de la TdP.

---

### 8. `TdP_fitter_4_1.py` — (Versión Histórica)

**Propósito**:  
Versión temprana del algoritmo de ajuste, desarrollada durante la fase exploratoria.

**Estado**:  
Archivado. No se usa en producción. Preserva la evolución conceptual del proyecto.

---

🔗 **Artículo asociado**: [Signatures-of-a-Fractal-Vacuum-in-Gaussian-Boson-Sampling](https://www.scribd.com/document/910269644/Signatures-of-a-Fractal-Vacuum-in-Gaussian-Boson-Sampling)

### `Model_Gaussian_Boson_Sampling.py` — Módulo de Adquisición y Preprocesamiento de Datos

Este script implementa el módulo central para el análisis de datos de muestreo bosónico gaussiano, diseñado para manejar los volúmenes masivos de datos esperados de experimentos como "Jiuzhang 4.0". Aunque conciso, es una obra maestra de ingeniería científica por su flexibilidad, eficiencia y robustez.

#### 1. **Flexibilidad y Robustez: `load_jiuzhang_data`**

La función `load_jiuzhang_data` está diseñada para ser **ag-nóstica al formato de entrada**, permitiendo la carga de datos desde:
- Archivos binarios (formato típico para grandes volúmenes de datos),
- Archivos CSV (para pruebas y validación),
- Arrays en memoria (para simulaciones como las realizadas en este proyecto).

Incluye manejo de errores con bloques `try-except` y excepciones `ValueError` para capturar problemas comunes (rutas incorrectas, formatos inválidos, dimensiones inconsistentes), garantizando que el análisis falle de forma **grácil y diagnóstica**, no de forma críptica.

#### 2. **Eficiencia Computacional: `compute_covariance`**

La función `compute_covariance` es clave para el rendimiento:

- **Matrices Dispersas (`scipy.sparse`)**: Para "Jiuzhang 4.0", la matriz de covarianza será de 8176 × 8176 (más de 66 millones de elementos). Al usar matrices dispersas, se evita el almacenamiento denso y se optimiza el cálculo.
  
- **Cálculo Vectorizado**: La fórmula:
  > (events_sparse.T @ events_sparse) / M - np.outer(mean_n, mean_n)
  
  es la forma más eficiente y numéricamente estable de calcular la matriz de covarianza, aprovechando operaciones de álgebra lineal optimizadas.

#### 3. **Modularidad y Claridad: `preprocess_jiuzhang`**

La función principal `preprocess_jiuzhang` encapsula todo el flujo de preprocesamiento:
- Toma datos crudos,
- Los valida,
- Calcula la matriz de eventos y la matriz de covarianza.

Es un **"caja negra" bien definida**, fácil de integrar en pipelines de análisis más grandes, como los necesarios para las tres pruebas de validación de la TdP.

#### 4. **Prueba Unitaria Incorporada**

El script incluye una **prueba unitaria automática** que:
- Genera datos sintéticos con propiedades conocidas,
- Ejecuta las funciones de carga y cálculo,
- Verifica que los resultados sean correctos.

Esta auto-validación es esencial: garantiza que el módulo funcione como se espera antes de aplicarlo a datos reales, brindando **confianza en la integridad del análisis**.

---

> Este módulo no solo procesa datos.  
> **Valida la conexión entre el experimento y la teoría.**  
> Es el primer eslabón en la cadena que podría confirmar que el universo tiene una estructura fractal p-ádica.
> 
## 🛠️ Requisitos

- Python 3.0
- `numpy`, `scipy`, `matplotlib`, `corner`
- `numba` (para aceleración JIT)
- `emcee` (para MCMC)
- `jupyter` (opcional, para notebooks)

```bash
pip install numpy scipy matplotlib numba emcee jupyter corner

