import numpy as np


def normalize(v: np.ndarray):
    norm = np.linalg.norm(v, axis=-1)
    if v.ndim == 2:
        return v / norm[:, None]
    return v / norm
