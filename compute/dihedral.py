import numpy as np
#from MDAnalysis.analysis.dihedrals import Dihedral

class Dihedral:
    def __init__(self, u, common="protein", atoms=[]):
        assert len(atoms)==4, "len(atoms) != 4"
        if isinstance(atoms[0], list):
            ag = u.select_atoms(common + " and resid %d and " %atoms[0][0] + "name %s" %atoms[0][1]) + \
                 u.select_atoms(common + " and resid %d and " %atoms[1][0] + "name %s" %atoms[1][1]) + \
                 u.select_atoms(common + " and resid %d and " %atoms[2][0] + "name %s" %atoms[2][1]) + \
                 u.select_atoms(common + " and resid %d and " %atoms[3][0] + "name %s" %atoms[3][1])

        else:
            ag = u.select_atoms(common + " and name %s" %atoms[0]) + \
                 u.select_atoms(common + " and name %s" %atoms[1]) + \
                 u.select_atoms(common + " and name %s" %atoms[2]) + \
                 u.select_atoms(common + " and name %s" %atoms[3])

        self.ag = ag
        self.u = u
        print("indices(0-based): ",  ag.indices)
        print("indices(1-based): ",  ag.indices+1)


    def run(self, b=0, e=None, skip=1, cossin=True):
        times = []
        angles = []
        for ts in self.u.trajectory[b:e:skip]:
            times.append(ts.time/1000)
            angles.append(self.ag.dihedral.value())

        angles = np.array(angles) * np.pi/180

        if cossin:
            output = np.array([times, np.cos(angles), np.sin(angles)])
        else:
            output = np.array([times, angles])

        return np.transpose(output)

            


