import numpy as np


class Camera:
    def __init__(self, dtype=np.float64):

        self.rotation_cg: np.ndarray = np.identity(3, dtype=dtype)
        self.translation_cg: np.ndarray = np.array([0, 0, 0], dtype=dtype)
