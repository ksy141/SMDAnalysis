from __future__ import absolute_import, division, print_function

import numpy as np
from MDAnalysis import Universe
from MDAnalysis import AtomGroup


class TrajectoryToNumpy:
    def __init__(self, u, selection='all'):
        assert isinstance(u, Universe)
        self.u     = u
        self.group = u.select_atoms(selection)
        assert isinstance(self.group, AtomGroup)

    @property
    def to_fac(self):
        """
        FAC order: FRAME/ATOM/COORD
        :return: FAC numpy array
        """
        fac = []
        for ts in self.u.trajectory:
            fac.append(self.group.positions)

        return np.array(fac)

    @property
    def to_afc(self):
        """
        AFC order: ATOM/FRAME/COORD
        :return: AFC numpy array
        """
        fac = self.to_fac
        afc = np.transpose(fac, (1, 0, 2))

        return afc
