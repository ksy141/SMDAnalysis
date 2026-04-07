from __future__ import absolute_import, division, print_function

import numpy as np

class Distance:
    def __init__(self, r1, r2, pbc):
        assert isinstance(r1, np.ndarray)
        assert isinstance(r2, np.ndarray)
        assert isinstance(pbc, np.ndarray)
        self.pbc = pbc
        self.r1  = r1
        self.r2  = r2

    def distance(self, pbc=None):
        dr = self.r1 - self.r2
        if pbc:
            dr -= self.pbc * np.round(dr/self.pbc)
            return np.linalg.norm(dr)
        else:
            return np.linalg.norm(dr)

    def distance2(self, pbc=None):
        dr = self.r1 - self.r2
        if pbc:
            dr -= self.pbc * np.round(dr/self.pbc)
            return np.sum(dr**2, axis=-1)
            #if dr.shape[1] == 3:
            #    return np.sum(dr**2, axis = 1)
            #else:
            #    return np.sum(dr**2)
        else:
            return np.sum(dr**2, axis=-1)
            #if dr.shape[1] == 3:
            #    return np.sum(dr**2, axis = 1)
            #else:
            #    return np.sum(dr**2)


