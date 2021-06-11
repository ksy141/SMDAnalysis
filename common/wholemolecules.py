from __future__ import absolute_import, division, print_function

import numpy as np
from MDAnalysis import AtomGroup

class Wholemolecules:
    '''
    Make wholemolecules of selection
    Wholemolecules(AtomGroup)
    '''
    def __init__(self, atomgroup, point=None):
        '''
        Make wholemolecules based on the first atom or point
        '''
        assert isinstance(atomgroup, AtomGroup)
        pbc = atomgroup.dimensions[0:3]
        positions = atomgroup.positions

        if point is None:
            dr  = positions - positions[0]
            dr -= pbc * np.round(dr/pbc)
            atomgroup.positions = positions[0] + dr

        else:
            dr  = positions - point
            dr -= pbc * np.round(dr/pbc)
            atomgroup.positions = point + dr




