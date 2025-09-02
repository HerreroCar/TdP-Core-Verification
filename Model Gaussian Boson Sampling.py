import numpy as np
from scipy import sparse
import pandas as pd

def load_jiuzhang_data(data_input, num_modes=144, format='binary'):
    """
    Carga datos Jiuzhang (CSV/binario o memoria). Retorna matriz eventos (M x N).
    Args:
        data_input: str (path a archivo) o np.ndarray (datos en memoria).
        num_modes: int, número de modos (e.g., 144 para Jiuzhang 3.0, 8176 para 4.0).
        format: str, 'binary' (bin file), 'csv' (CSV file), o 'memory' (np.ndarray).
    Returns:
        np.ndarray, matriz eventos (M x N).
    Raises:
        ValueError: Si formato o input son inválidos.
    """
    if format == 'memory':
        if not isinstance(data_input, np.ndarray):
            raise ValueError("Para format='memory', data_input debe ser np.ndarray")
        data = data_input
        if data.shape[1] != num_modes:
            raise ValueError(f"Data shape[1]={data.shape[1]} no coincide con num_modes={num_modes}")
    elif format == 'binary':
        if not isinstance(data_input, str):
            raise ValueError("Para format='binary', data_input debe ser un path str")
        try:
            data = np.fromfile(data_input, dtype=np.uint8).reshape(-1, num_modes)
        except ValueError as e:
            raise ValueError(f"Error al leer binario: {e}. Verifica num_modes o archivo.")
    elif format == 'csv':
        if not isinstance(data_input, str):
            raise ValueError("Para format='csv', data_input debe ser un path str")
        try:
            data = pd.read_csv(data_input, header=None).values.astype(np.uint8)
        except Exception as e:
            raise ValueError(f"Error al leer CSV: {e}. Verifica formato del archivo.")
    else:
        raise ValueError(f"Formato no soportado: {format}")
    return data

def compute_covariance(events):
    """
    Calcula matriz covarianza Σ_ij = <n_i n_j> - <n_i><n_j>.
    Usa sparse para N grande (e.g., 8176).
    Args:
        events: np.ndarray, matriz M x N (eventos x modos).
    Returns:
        scipy.sparse.csr_matrix, matriz covarianza N x N.
    """
    M = events.shape[0]
    mean_n = np.mean(events, axis=0)
    # Usa sparse matrix para eficiencia en N grande
    events_sparse = sparse.csr_matrix(events)
    cov = (events_sparse.T @ events_sparse) / M - np.outer(mean_n, mean_n)
    return sparse.csr_matrix(cov)

def preprocess_jiuzhang(data_input, num_modes=144, format='binary'):
    """
    Preprocesa datos Jiuzhang. Retorna eventos y covarianza.
    Args:
        data_input: str o np.ndarray, datos crudos o en memoria.
        num_modes: int, número de modos.
        format: str, 'binary', 'csv', o 'memory'.
    Returns:
        tuple, (eventos: np.ndarray, covarianza: scipy.sparse.csr_matrix).
    """
    events = load_jiuzhang_data(data_input, num_modes, format)
    cov = compute_covariance(events)
    return events, cov

# Test con datos sintéticos (mimic Jiuzhang 3.0, 144 modos, 10k eventos)
np.random.seed(42)
synthetic_data = np.random.poisson(lam=2, size=(10000, 144)).astype(np.uint8)
events, cov = preprocess_jiuzhang(synthetic_data, num_modes=144, format='memory')
print("Eventos shape:", events.shape)
print("Covarianza shape:", cov.shape)